
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <util/misc/formio.h>
#include <util/group/memory.h>
#include <util/group/memory.cc>

#include <util/group/memproc.h>

#ifdef HAVE_NX
#  ifdef HAVE_HRECV
#    include <util/group/mempgon.h>
#  endif
#  include <util/group/messpgon.h>
#  include <util/group/memipgon.h>
#endif

#ifdef HAVE_SYSV_IPC
#  include <util/group/messshm.h>
#  include <util/group/memshm.h>
#endif

#if defined(HAVE_MPL)
#  include <util/group/memmpl.h>
#endif
#if defined(HAVE_MPI)
#  include <util/group/messmpi.h>
#  include <util/group/memmpi.h>
#endif

//////////////////////////////////////////////////////////////////////
// MemoryGrpBuf template instantiations

#ifdef __GNUG__
template class MemoryGrpBuf<double>;
template class MemoryGrpBuf<int>;
template class MemoryGrpBuf<char>;
template class MemoryGrpBuf<unsigned char>;
#endif

//////////////////////////////////////////////////////////////////////
// MemoryGrp members

#define CLASSNAME MemoryGrp
#define PARENTS public DescribedClass
#include <util/class/classia.h>
void *
MemoryGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

MemoryGrp::MemoryGrp()
{
  use_locks_ = 1;
  debug_ = 0;
}

MemoryGrp::MemoryGrp(const RefKeyVal& keyval)
{
  use_locks_ = 1;
  debug_ = keyval->intvalue("debug");
}

MemoryGrp::~MemoryGrp()
{
  delete[] offsets_;
}

MemoryGrp *
MemoryGrp::initial_memorygrp()
{
  int argc = 0;
  initial_memorygrp(argc,0);
}

MemoryGrp *
MemoryGrp::initial_memorygrp(int &argc, char *argv[])
{
  MemoryGrp *grp = 0;

  char *keyval_string = 0;

  // see if a memory group is given on the command line
  if (argc && argv) {
      for (int i=0; i<argc; i++) {
	  if (argv[i] && !strcmp(argv[i], "-memorygrp")) {
              i++;
              if (i >= argc) {
                  cerr << "-memorygrp must be following by an argument"
                       << endl;
                  abort();
                }
              keyval_string = argv[i];
              // permute the memorygrp arguments to the end of argv
              char *tmp = argv[argc-2];
              argv[argc-2] = argv[i-1];
              argv[i-1] = tmp;
              tmp = argv[argc-1];
              argv[argc-1] = argv[i];
              argv[i] = tmp;
              break;
            }
        }
    }

  if (!keyval_string) {
      // find out if the environment gives the containing memory group
      keyval_string = getenv("MEMORYGRP");
      if (keyval_string) {
          if (!strncmp("MEMORYGRP=", keyval_string, 11)) {
              keyval_string = strchr(keyval_string, '=');
            }
          if (*keyval_string == '=') keyval_string++;
        }
    }

  // if keyval input for a memory group was found, then
  // create it.
  if (keyval_string) {
      //cout << "Creating MemoryGrp from \"" << keyval_string << "\"" << endl;
      RefParsedKeyVal strkv = new ParsedKeyVal();
      strkv->parse_string(keyval_string);
      RefDescribedClass dc = strkv->describedclassvalue();
      grp = MemoryGrp::castdown(dc.pointer());
      if (dc.null()) {
          cerr << "initial_memorygrp: couldn't find a MemoryGrp in "
               << keyval_string << endl;
          abort();
        }
      else if (!grp) {
          cerr << "initial_memorygrp: wanted MemoryGrp but got "
               << dc->class_name() << endl;
          abort();
        }
      // prevent an accidental delete
      grp->reference();
      strkv = 0;
      dc = 0;
      // accidental delete not a problem anymore since all smart pointers
      // to grp are dead
      grp->dereference();
      return grp;
    }

  RefMessageGrp msg = MessageGrp::get_default_messagegrp();
  if (msg.null()) {
      cerr << scprintf("MemoryGrp::create_memorygrp: requires default msg\n");
      abort();
    }
#ifdef HAVE_NX
  else if (msg->class_desc() == ParagonMessageGrp::static_class_desc()) {
#ifdef HAVE_HRECV
      grp = new ParagonMemoryGrp(msg);
#else
      grp = new IParagonMemoryGrp(msg);
#endif
    }
#endif
#if defined(HAVE_MPL)
  else if (msg->class_desc() == MPIMessageGrp::static_class_desc()) {
      grp = new MPLMemoryGrp(msg);
    }
#endif
#if defined(HAVE_MPI)
  else if (msg->class_desc() == MPIMessageGrp::static_class_desc()) {
      grp = new MPIMemoryGrp(msg);
    }
#endif
#ifdef HAVE_SYSV_IPC
  else if (msg->class_desc() == ShmMessageGrp::static_class_desc()) {
      grp = new ShmMemoryGrp(msg);
    }
#endif
  else if (msg->n() == 1) {
      grp = new ProcMemoryGrp();
    }
  else {
      cerr << scprintf("MemoryGrp::create_memorygrp: cannot create "
              "default for \"%s\"\n.", msg->class_name());
      abort();
    }

  if (!grp) {
      cerr << scprintf("WARNING: MemoryGrp::initial_memorygrp(): failed\n");
    }

  return grp;
}

void
MemoryGrp::lock(int b)
{
  use_locks_ = b;
}

void
MemoryGrp::release_read_(int offset, int size)
{
  if (offset < offsets_[me_]
      || offset + size > offsets_[me_+1]) {
      cerr << scprintf("MemoryGrp::release_read_: bad args\n");
      abort();
    }

  locks_.decrement(offset, offset + size);
}

void
MemoryGrp::release_write_(int offset, int size)
{
  if (offset < offsets_[me_]
      || offset + size > offsets_[me_+1]) {
      cerr << scprintf("MemoryGrp::release_write_: bad args\n");
      abort();
    }

  locks_.increment(offset, offset + size);
}

// return 1 if the lock was obtained, otherwise 0
int
MemoryGrp::obtain_read_(int offset, int size)
{
  if (offset < offsets_[me_]
      || offset + size > offsets_[me_+1]) {
      cerr << scprintf("MemoryGrp::obtain_read_: bad args\n");
      abort();
    }

  if (!locks_.checkgr(offset, offset + size, -1)) return 0;
  locks_.increment(offset, offset + size);
  return 1;
}

// return 1 if the lock was obtained, otherwise 0
int
MemoryGrp::obtain_write_(int offset, int size)
{
  if (offset < offsets_[me_]
      || offset + size > offsets_[me_+1]) {
      cerr << scprintf("MemoryGrp::obtain_write_: bad args\n");
      abort();
    }

  if (!locks_.checkeq(offset, offset + size, 0)) return 0;
  locks_.decrement(offset, offset + size);
  return 1;
}

void
MemoryGrp::activate()
{
}

void
MemoryGrp::deactivate()
{
}

void
MemoryGrp::print(ostream&o)
{
  o << scprintf("MemoryGrp (node %d):\n", me());
  locks_.print(o);
  o << scprintf("%d: n = %d\n", me(), n());
  for (int i=0; i<=n_; i++) {
      o << scprintf("%d: offset[%d] = %5d\n", me(), i, offsets_[i]);
    }
}

void *
MemoryGrp::obtain_writeonly(int offset, int size)
{
  return obtain_readwrite(offset, size);
}

void
MemoryGrp::sum_reduction(double *data, int doffset, int dlength)
{
  int offset = doffset * sizeof(double);
  int length = dlength * sizeof(double);

  if (offset + length > totalsize()) {
      cerr << scprintf("MemoryGrp::sum_reduction: arg out of range\n");
      abort();
    }

  double *source_data = (double*) obtain_readwrite(offset, length);

  for (int i=0; i<dlength; i++) {
      source_data[i] += data[i];
    }

  release_write((void*) source_data, offset, length);
}

void
MemoryGrp::sum_reduction_on_node(double *data, int doffset, int dlength,
                                 int node)
{
  if (node == -1) node = me();

  sum_reduction(data, doffset + offset(node)/sizeof(double),
                dlength);
}
