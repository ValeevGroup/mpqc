
#ifndef _util_group_memmsg_cc
#define _util_group_memmsg_cc

#ifdef __GNUC__
#pragma implementation
#endif

#include <util/group/memmsg.h>

#define CLASSNAME MsgMemoryGrp
#define PARENTS public MemoryGrp
#include <util/class/classia.h>
void *
MsgMemoryGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  MemoryGrp::_castdown(cd);
  return do_castdowns(casts,cd);
}

MsgMemoryGrp::MsgMemoryGrp(const RefMessageGrp &msg)
{
  msg_ = msg;
  n_ = msg->n();
  me_ = msg->me();

  offsets_ = 0;
}

MsgMemoryGrp::MsgMemoryGrp(const RefKeyVal &keyval):
  MemoryGrp(keyval)
{
  RefMessageGrp msg = keyval->describedclassvalue("message");
  if (msg.null()) {
      msg = MessageGrp::get_default_messagegrp();
    }
  if (msg.null()) {
      cerr << "MsgMemoryGrp(const RefKeyVal&): couldn't find MessageGrp"
           << endl;
      abort();
    }

  msg_ = msg;
  n_ = msg->n();
  me_ = msg->me();

  offsets_ = 0;
}

MsgMemoryGrp::~MsgMemoryGrp()
{
  delete[] offsets_;
}

void
MsgMemoryGrp::set_localsize(int localsize)
{
  delete[] offsets_;

  offsets_ = new int[n_ + 1];
  int *sizes = new int[n_];

  int i;
  for (i=0; i<n_; i++) sizes[i] = 0;
  sizes[me_] = localsize;

  msg_->sum(sizes, n_);

  offsets_[0] = 0;
  for (i=1; i<=n_; i++) {
      offsets_[i] = sizes[i-1] + offsets_[i-1];
    }

  delete[] sizes;
}

void
MsgMemoryGrp::sync()
{
  msg_->sync();
}

#endif
