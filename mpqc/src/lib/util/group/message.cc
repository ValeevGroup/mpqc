
#ifdef __GNUC__
#pragma implementation
#endif

#include <stdio.h>
#include <string.h>

#include <util/group/message.h>
#include <util/group/topology.h>
#include <util/group/hcube.h>
#include <util/class/classMap.h>

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

MessageGrp::MessageGrp(const RefKeyVal&):
  me_(-1),
  n_(-1),
  classdesc_to_index_(-1),
  index_to_classdesc_(0)
{
}

MessageGrp::MessageGrp():
  me_(-1),
  n_(-1),
  classdesc_to_index_(-1),
  index_to_classdesc_(0)
{
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
                  printf("node %d adding class %d = \"%s\"\n",
                         me, iclass, classes.contents(J)->name());
#endif
                  strcpy(currentbuffer,classes.contents(J)->name());
                  currentbuffer += strlen(classes.contents(J)->name()) + 1;
                }
            }
#ifdef DEBUG
          printf("node %d bcast n_new_class = %d\n",me,n_new_class);
#endif
          bcast(&n_new_class,1,i);
#ifdef DEBUG
          printf("node %d finished bcast\n",me);
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
          printf("node %d begin recv bcast\n",me);
#endif
          bcast(&n_new_class,1,i);
#ifdef DEBUG
          printf("node %d recv bcast n_new_class = %d\n",me,n_new_class);
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
                      printf("node %d adding \"%s\"\n", me, currentbuffer);
#endif
                      classdesc_to_index_[cd] = iclass;
                    }
#ifdef DEBUG
                  else {
                      printf("node %d failed to add \"%s\"\n",
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
      printf("MessageGrp: registered %d classes and %d nodes\n",
             nclass_, n_);
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
  RefGlobalMsgIter i(topology_->global_msg_iter(this, from));
  for (i->forwards(); !i->done(); i->next()) {
      if (i->send()) {
          raw_send(i->sendto(), data, nbyte);
        }
      if (i->recv()) {
          raw_recv(i->recvfrom(), data, nbyte);
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
