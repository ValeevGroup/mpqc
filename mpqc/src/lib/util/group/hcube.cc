//
// hcube.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
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

#include <util/group/message.h>
#include <util/group/topology.h>
#include <util/group/hcube.h>

using namespace sc;

static ClassDesc HypercubeGMI_cd(
  typeid(HypercubeGMI),"HypercubeGMI",1,"public GlobalMsgIter",
  0, 0, 0);

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

static ClassDesc HypercubeTopology_cd(
  typeid(HypercubeTopology),"HypercubeTopology",1,"public MachineTopology",
  0, create<HypercubeTopology>, 0);

HypercubeTopology::HypercubeTopology()
{
}

HypercubeTopology::HypercubeTopology(const Ref<KeyVal>& keyval):
  MachineTopology(keyval)
{
}

HypercubeTopology::~HypercubeTopology()
{
}

Ref<GlobalMsgIter>
HypercubeTopology::global_msg_iter(const Ref<MessageGrp>& grp,
                                   int root)
{
  return new HypercubeGMI(grp->n(), grp->me(), root);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
