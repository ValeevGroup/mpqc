
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <util/group/memory.h>
#include <util/group/memory.cc>

#include <util/group/memproc.h>

#ifdef HAVE_HRECV
#  include <util/group/messpgon.h>
#  include <util/group/mempgon.h>
#endif

#ifdef HAVE_SYSV_IPC
#  include <util/group/messshm.h>
#  include <util/group/memshm.h>
#endif

#if defined(HAVE_MPL) && defined(HAVE_MPI)
#  include <util/group/memmpl.h>
#  include <util/group/messmpi.h>
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
}

MemoryGrp::~MemoryGrp()
{
  delete[] offsets_;
}

MemoryGrp *
MemoryGrp::create_memorygrp(int localsize)
{
  MemoryGrp *ret = 0;
  RefMessageGrp msg = MessageGrp::get_default_messagegrp();
  if (msg.null()) {
      fprintf(stderr, "MemoryGrp::create_memorygrp: requires default msg\n");
      abort();
    }
#ifdef HAVE_HRECV
  else if (msg->class_desc() == ParagonMessageGrp::static_class_desc()) {
      ret = new ParagonMemoryGrp(msg, localsize);
    }
#endif
#if defined(HAVE_MPL) && defined(HAVE_MPI)
  else if (msg->class_desc() == MPIMessageGrp::static_class_desc()) {
      printf("creating mplmemorygrp\n");
      ret = new MPLMemoryGrp(msg, localsize);
    }
#endif
#ifdef HAVE_SYSV_IPC
  else if (msg->class_desc() == ShmMessageGrp::static_class_desc()) {
      ret = new ShmMemoryGrp(msg, localsize);
    }
#endif
  else if (msg->n() == 1) {
      ret = new ProcMemoryGrp(localsize);
    }
  else {
      fprintf(stderr, "MemoryGrp::create_memorygrp: cannot create "
              "default for \"%s\"\n.", msg->class_name());
      abort();
    }

  if (!ret) {
      fprintf(stderr, "WARNING: MemoryGrp::create_memorygrp(): failed\n");
    }

  return ret;
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
      fprintf(stderr,"MemoryGrp::release_read_: bad args\n");
      abort();
    }

  locks_.decrement(offset, offset + size);
}

void
MemoryGrp::release_write_(int offset, int size)
{
  if (offset < offsets_[me_]
      || offset + size > offsets_[me_+1]) {
      fprintf(stderr,"MemoryGrp::release_write_: bad args\n");
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
      fprintf(stderr,"MemoryGrp::obtain_read_: bad args\n");
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
      fprintf(stderr,"MemoryGrp::obtain_write_: bad args\n");
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
MemoryGrp::print(FILE *fp)
{
  fprintf(fp, "MemoryGrp (node %d):\n", me());
  locks_.print(fp);
  fprintf(fp, "%d: n = %d\n", me(), n());
  for (int i=0; i<=n_; i++) {
      fprintf(fp, "%d: offset[%d] = %5d\n", me(), i, offsets_[i]);
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
      fprintf(stderr, "MemoryGrp::sum_reduction: arg out of range\n");
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
