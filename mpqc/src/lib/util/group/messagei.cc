//
// messagei.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <string.h>

#include <util/misc/formio.h>

#include <util/group/message.h>

#include <util/group/topology.h>
#include <util/group/hcube.h>
#include <util/class/classMap.h>
#ifdef HAVE_NX
#  include <util/group/messpgon.h>
#endif
#ifdef HAVE_MPI
#  include <util/group/messmpi.h>
#endif

#define CLASSNAME MessageGrp
#define PARENTS public DescribedClass
#include <util/class/classia.h>
void *
MessageGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

MessageGrp::MessageGrp(const RefKeyVal& keyval):
  me_(-1),
  n_(-1),
  classdesc_to_index_(-1),
  index_to_classdesc_(0)
{
  gop_max_ = keyval->intvalue("gop_max");
  if (keyval->error() != KeyVal::OK) gop_max_ = 320000;
  debug_ = keyval->booleanvalue("debug");
}

MessageGrp::MessageGrp():
  me_(-1),
  n_(-1),
  classdesc_to_index_(-1),
  index_to_classdesc_(0)
{
  gop_max_ = 320000;
  debug_ = 0;
}

MessageGrp::~MessageGrp()
{
  if (index_to_classdesc_) delete[] index_to_classdesc_;
}

static RefMessageGrp default_messagegrp;

void
MessageGrp::set_default_messagegrp(const RefMessageGrp& grp)
{
  default_messagegrp = grp;
}

MessageGrp*
MessageGrp::get_default_messagegrp()
{
  if (default_messagegrp.null()) {
      default_messagegrp = new ProcMessageGrp;
    }
  return default_messagegrp.pointer();
}

MessageGrp *
MessageGrp::initial_messagegrp(int &argc, char** argv)
{
  MessageGrp *grp = 0;

  char *keyval_string = 0;

  // see if a message group is given on the command line
  if (argc && argv) {
      for (int i=0; i<argc; i++) {
	  if (argv[i] && !strcmp(argv[i], "-messagegrp")) {
              i++;
              if (i >= argc) {
                  cerr << "-messagegrp must be following by an argument"
                       << endl;
                  abort();
                }
              keyval_string = argv[i];
              // permute the messagegrp arguments to the end of argv
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
      // find out if the environment gives the containing message group
      keyval_string = getenv("MESSAGEGRP");
      if (keyval_string) {
          if (!strncmp("MESSAGEGRP=", keyval_string, 11)) {
              keyval_string = strchr(keyval_string, '=');
            }
          if (*keyval_string == '=') keyval_string++;
        }
    }

  // if keyval input for a message group was found, then
  // create it.
  if (keyval_string) {
      if (keyval_string[0] == '\0') return 0;
      //cout << "Creating MessageGrp from \"" << keyval_string << "\"" << endl;
      RefParsedKeyVal strkv = new ParsedKeyVal();
      strkv->parse_string(keyval_string);
      RefDescribedClass dc = strkv->describedclassvalue();
      grp = MessageGrp::castdown(dc.pointer());
      if (dc.null()) {
          cerr << "initial_messagegrp: couldn't find a MessageGrp in "
               << keyval_string << endl;
          abort();
        }
      else if (!grp) {
          cerr << "initial_messagegrp: wanted MessageGrp but got "
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

  // if certain libraries have been compiled in, use those message groups
#if defined(HAVE_PUMA_MPI2)
  grp = new MPIMessageGrp;
  if (grp->n() == 1) { delete grp; return new ProcMessageGrp; }
  else return grp;
#elif defined(HAVE_NX)
  grp = new ParagonMessageGrp;
  if (grp->n() == 1) { delete grp; return new ProcMessageGrp; }
  else return grp;
#endif
  return 0;
}

void
MessageGrp::initialize(int me, int n)
{
  // This member is called by a CTOR and, ultimately, causes
  // 'this' to be converted into a temporary Ref which causes
  // this to be deleted (very bad), so reference 'this'
  // (and dereference this down below).
  this->reference();

  if (topology_.null()) {
      topology_ = new HypercubeTopology();
    }

  int i;
  Pix J;
  
  me_ = me;
  n_ = n;

  // get all of the classes known on this node
  ClassKeyClassDescPMap& classes = ClassDesc::all();

  // Keeps count of how many classes are known.
  int iclass = 0;

  for (i=0; i<n; i++) {
      if (i==me) {
          // Find out how many of my classes are not yet in the
          // classdesc to index map.
          int n_new_class = 0;
          int buffer_size = 0;
          for (J=classes.first(); J; classes.next(J)) {
              // I'm not sure where the null entries come from, but
              // they are bad--skip them.
              if (classes.contents(J) == 0) continue;
              if (!classdesc_to_index_.contains(classes.contents(J))) {
                  n_new_class++;
                  buffer_size += strlen(classes.contents(J)->name()) + 1;
                }
            }
          char* buffer = new char[buffer_size];
          char* currentbuffer = buffer;
          for (J=classes.first(); J; classes.next(J)) {
              if (classes.contents(J) == 0) continue;
              if (!classdesc_to_index_.contains(classes.contents(J))) {
                  classdesc_to_index_[classes.contents(J)] = iclass;
                  iclass++;
#ifdef DEBUG
                  cout << scprintf("node %d adding class %d = \"%s\"\n",
                                   me, iclass, classes.contents(J)->name());
#endif
                  strcpy(currentbuffer,classes.contents(J)->name());
                  currentbuffer += strlen(classes.contents(J)->name()) + 1;
                }
            }
#ifdef DEBUG
          cout << scprintf("node %d bcast n_new_class = %d\n",me,n_new_class);
#endif
          bcast(&n_new_class,1,i);
#ifdef DEBUG
          cout << scprintf("node %d finished bcast\n",me);
#endif
          if (n_new_class) {
              bcast(&buffer_size,1,i);
              bcast(buffer,buffer_size,i);
            }
          delete[] buffer;
        }
      else {
          int j;
          // Get new classnames and indices from node i.
          int n_new_class;
#ifdef DEBUG
          cout << scprintf("node %d begin recv bcast\n",me);
#endif
          bcast(&n_new_class,1,i);
#ifdef DEBUG
          cout << scprintf("node %d recv bcast n_new_class = %d\n",
                           me,n_new_class);
#endif
          if (n_new_class) {
              int buffer_size;
              bcast(&buffer_size,1,i);
              char* buffer = new char[buffer_size];
              char* currentbuffer = buffer;
              bcast(buffer,buffer_size,i);
              for (j=0; j<n_new_class; j++) {
                  ClassDescP cd = ClassDesc::name_to_class_desc(currentbuffer);
                  if (cd) {
#ifdef DEBUG
                      cout << scprintf("node %d adding \"%s\"\n",
                                       me, currentbuffer);
#endif
                      classdesc_to_index_[cd] = iclass;
                    }
#ifdef DEBUG
                  else {
                      cout << scprintf("node %d failed to add \"%s\"\n",
                                       me, currentbuffer);
                    }
#endif
                  iclass++;
                  // advance the currentbuffer to the next classname
                  while(*currentbuffer != '\0') currentbuffer++;
                  currentbuffer++;
                }
              delete[] buffer;
            }
        }
    }
  nclass_ = iclass;

  // Construct the mapping of index to classdesc.
  index_to_classdesc_ = new ClassDescP[iclass];
  for (i=0; i<nclass_; i++) {
      index_to_classdesc_[i] = 0;
    }
  for (J=classes.first(); J; classes.next(J)) {
      if (classdesc_to_index_.contains(classes.contents(J))) {
          index_to_classdesc_[classdesc_to_index_[classes.contents(J)]] = classes.contents(J);
        }
    }

  if (me_ == 0) {
    cout << "MessageGrp(" << class_name() << "): registered "
         << nclass_ << " classes and " << n_ << " nodes" << endl;
    }

  this->dereference();
}

// Sequential send routines

void
MessageGrp::send(int target, double* data, int ndata)
{
  raw_send(target, data, ndata*sizeof(double));
}
void
MessageGrp::send(int target, short* data, int ndata)
{
  raw_send(target, data, ndata*sizeof(short));
}
void
MessageGrp::send(int target, long* data, int ndata)
{
  raw_send(target, data, ndata*sizeof(long));
}
void
MessageGrp::send(int target, float* data, int ndata)
{
  raw_send(target, data, ndata*sizeof(float));
}
void
MessageGrp::send(int target, int* data, int ndata)
{
  raw_send(target, data, ndata*sizeof(int));
}
void
MessageGrp::send(int target, char* data, int ndata)
{
  raw_send(target, data, ndata);
}
void
MessageGrp::send(int target, unsigned char* data, int ndata)
{
  raw_send(target, data, ndata);
}
void
MessageGrp::send(int target, signed char* data, int ndata)
{
  raw_send(target, data, ndata);
}

// Sequential receive routines

void
MessageGrp::recv(int sender, double* data, int ndata)
{
  raw_recv(sender, data, ndata*sizeof(double));
}
void
MessageGrp::recv(int sender, short* data, int ndata)
{
  raw_recv(sender, data, ndata*sizeof(short));
}
void
MessageGrp::recv(int sender, long* data, int ndata)
{
  raw_recv(sender, data, ndata*sizeof(long));
}
void
MessageGrp::recv(int sender, float* data, int ndata)
{
  raw_recv(sender, data, ndata*sizeof(float));
}
void
MessageGrp::recv(int sender, int* data, int ndata)
{
  raw_recv(sender, data, ndata*sizeof(int));
}
void
MessageGrp::recv(int sender, char* data, int ndata)
{
  raw_recv(sender, data, ndata);
}
void
MessageGrp::recv(int sender, unsigned char* data, int ndata)
{
  raw_recv(sender, data, ndata);
}
void
MessageGrp::recv(int sender, signed char* data, int ndata)
{
  raw_recv(sender, data, ndata);
}

// Typed send routines

void
MessageGrp::sendt(int target, int type, double* data, int ndata)
{
  raw_sendt(target, type, data, ndata*sizeof(double));
}
void
MessageGrp::sendt(int target, int type, short* data, int ndata)
{
  raw_sendt(target, type, data, ndata*sizeof(short));
}
void
MessageGrp::sendt(int target, int type, long* data, int ndata)
{
  raw_sendt(target, type, data, ndata*sizeof(long));
}
void
MessageGrp::sendt(int target, int type, float* data, int ndata)
{
  raw_sendt(target, type, data, ndata*sizeof(float));
}
void
MessageGrp::sendt(int target, int type, int* data, int ndata)
{
  raw_sendt(target, type, data, ndata*sizeof(int));
}
void
MessageGrp::sendt(int target, int type, char* data, int ndata)
{
  raw_sendt(target, type, data, ndata);
}
void
MessageGrp::sendt(int target, int type, unsigned char* data, int ndata)
{
  raw_sendt(target, type, data, ndata);
}
void
MessageGrp::sendt(int target, int type, signed char* data, int ndata)
{
  raw_sendt(target, type, data, ndata);
}

// Typed receive routines

void
MessageGrp::recvt(int type, double* data, int ndata)
{
  raw_recvt(type, data, ndata*sizeof(double));
}
void
MessageGrp::recvt(int type, short* data, int ndata)
{
  raw_recvt(type, data, ndata*sizeof(short));
}
void
MessageGrp::recvt(int type, long* data, int ndata)
{
  raw_recvt(type, data, ndata*sizeof(long));
}
void
MessageGrp::recvt(int type, float* data, int ndata)
{
  raw_recvt(type, data, ndata*sizeof(float));
}
void
MessageGrp::recvt(int type, int* data, int ndata)
{
  raw_recvt(type, data, ndata*sizeof(int));
}
void
MessageGrp::recvt(int type, char* data, int ndata)
{
  raw_recvt(type, data, ndata);
}
void
MessageGrp::recvt(int type, unsigned char* data, int ndata)
{
  raw_recvt(type, data, ndata);
}
void
MessageGrp::recvt(int type, signed char* data, int ndata)
{
  raw_recvt(type, data, ndata);
}

// Broadcast operations

void
MessageGrp::bcast(double*data, int ndata, int from)
{
  raw_bcast(data, ndata*sizeof(double), from);
}
void
MessageGrp::bcast(short*data, int ndata, int from)
{
  raw_bcast(data, ndata*sizeof(short), from);
}
void
MessageGrp::bcast(long*data, int ndata, int from)
{
  raw_bcast(data, ndata*sizeof(long), from);
}
void
MessageGrp::bcast(float*data, int ndata, int from)
{
  raw_bcast(data, ndata*sizeof(float), from);
}
void
MessageGrp::bcast(int*data, int ndata, int from)
{
  raw_bcast(data, ndata*sizeof(int), from);
}
void
MessageGrp::bcast(char*data, int ndata, int from)
{
  raw_bcast(data, ndata, from);
}
void
MessageGrp::bcast(unsigned char*data, int ndata, int from)
{
  raw_bcast(data, ndata, from);
}
void
MessageGrp::bcast(signed char*data, int ndata, int from)
{
  raw_bcast(data, ndata, from);
}

// Global classdesc indices

int
MessageGrp::classdesc_to_index(const ClassDesc* cdptr)
{
  if (classdesc_to_index_.contains((ClassDesc*)cdptr)) {
      return classdesc_to_index_[(ClassDesc*)cdptr];
    }
  else {
      return -1;
    }
}

const ClassDesc*
MessageGrp::index_to_classdesc(int index)
{
  if (index < 0 || index >= nclass_) {
      return 0;
    }
  else {
      return index_to_classdesc_[index];
    }
}

void
MessageGrp::raw_bcast(void* data, int nbyte, int from)
{
  int nbyte_actual = nbyte;
  int tgop_max = nbyte;
  if (gop_max_ != 0) {
      tgop_max = gop_max_;
      gop_max_ = 0;
      bcast(nbyte_actual,from);
      gop_max_ = tgop_max;
    }
  for (int idat=0; idat<nbyte_actual; idat+=tgop_max) {
      int ndat = (idat+tgop_max>nbyte_actual)?(nbyte_actual-idat):tgop_max;
      RefGlobalMsgIter i(topology_->global_msg_iter(this, from));
      for (i->forwards(); !i->done(); i->next()) {
          if (i->send()) {
              raw_send(i->sendto(), &((char*)data)[idat], ndat);
            }
          if (i->recv()) {
              raw_recv(i->recvfrom(), &((char*)data)[idat], ndat);
            }
        }
    }
}

void
MessageGrp::sync()
{
  RefGlobalMsgIter i(topology_->global_msg_iter(this, 0));

  for (i->backwards(); !i->done(); i->next()) {
      if (i->send()) {
          raw_send(i->sendto(), 0, 0);
        }
      if (i->recv()) {
          raw_recv(i->recvfrom(), 0, 0);
        }
    }
  
  for (i->forwards(); !i->done(); i->next()) {
      if (i->send()) {
          raw_send(i->sendto(), 0, 0);
        }
      if (i->recv()) {
          raw_recv(i->recvfrom(), 0, 0);
        }
    }
}

void
MessageGrp::collect(const double *part, const int *lengths, double *whole)
{
  raw_collect(part,lengths,whole,sizeof(double));
}

void
MessageGrp::raw_collect(const void *part, const int *lengths, void *whole,
                        int bytes_per_datum)
{
  int offset = 0;
  for (int i=0; i<n_; i++) {
      int nbytes = lengths[i];
      if (i==me_) memcpy(&((char*)whole)[offset], part, nbytes);
      raw_bcast(&((char*)whole)[offset], nbytes, i);
      offset += nbytes;
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
