
#include "message.h"
#include "topology.h"
#include "hcube.h"

#define CLASSNAME HypercubeGMI
#define PARENTS public GlobalMsgIter
#include <util/class/classi.h>
void *
HypercubeGMI::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = GlobalMsgIter::_castdown(cd);
  return do_castdowns(casts,cd);
}

HypercubeGMI::HypercubeGMI(int nproc, int me, int root):
  GlobalMsgIter(nproc, me, root)
{
  int i;

  // compute the number of steps needed
  i = nproc;
  n_ = 0;
  while(i>1) {
      i = i>>1;
      n_++;
    }
  if (1<<n_ != nproc) n_++;
  if (root != 0) n_++;

  // compute the size of the hypercube that embeds these processors
  nhyper_ = 1;
  while (nhyper_ < nproc) nhyper_ = nhyper_ << 1;
}

HypercubeGMI::~HypercubeGMI()
{
}

int
HypercubeGMI::fwdsendto()
{
  int offset;
  if (root_ != 0) {
      if (i_ == 0) {
          if (root_ == me_) return 0;
          return -1;
        }
      else offset = 1;
    }
  else offset = 0;

  int bit = 1<<(i_-offset);
  int highbits = (nhyper_-2)<<(i_-offset);
  // if i don't have this bit and none above, then i'm a sender
  if (!(me_&bit) && !(me_&highbits)) {
      int target = me_ + bit;
      if (target >= nproc_) return -1;
      // already got this one
      if (target == root_) return -1;
      return target;
    }
  return -1;
}

int
HypercubeGMI::fwdsend()
{
  return fwdsendto() != -1;
}

int
HypercubeGMI::fwdrecvfrom()
{
  int offset;
  if (root_ != 0) {
      if (i_ == 0) {
          if (me_ == 0) return root_;
          return -1;
        }
      else offset = 1;
    }
  else offset = 0;

  // already got this one
  if (me_ == root_) return -1;

  int bit = 1<<(i_-offset);
  int highbits = (nhyper_-2)<<(i_-offset);
  if (!(me_&bit) && !(me_&highbits)) {
      return -1;
    }
  else if (!(me_&highbits)) {
      int source = me_ - bit;
      return source;
    }
  return -1;
}

int
HypercubeGMI::fwdrecv()
{
  return fwdrecvfrom() != -1;
}

///////////////////////////////////////////////////////////////////////////
// HypercubeTopology members

#define CLASSNAME HypercubeTopology
#define PARENTS public MachineTopology
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
HypercubeTopology::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = MachineTopology::_castdown(cd);
  return do_castdowns(casts,cd);
}

HypercubeTopology::HypercubeTopology()
{
}

HypercubeTopology::HypercubeTopology(const RefKeyVal& keyval):
  MachineTopology(keyval)
{
}

HypercubeTopology::~HypercubeTopology()
{
}

RefGlobalMsgIter
HypercubeTopology::global_msg_iter(const RefMessageGrp& grp,
                                   int root)
{
  return new HypercubeGMI(grp->n(), grp->me(), root);
}
