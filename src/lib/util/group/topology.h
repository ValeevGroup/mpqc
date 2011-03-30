//
// topology.h
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

#ifndef _util_group_topology_h
#define _util_group_topology_h

#include <util/class/class.h>
#include <util/keyval/keyval.h>

namespace sc {

class GlobalMsgIter: public DescribedClass {
  protected:
    int me_;
    int nproc_;
    int root_;
    int i_;
    int n_; // the number of steps--intialized by derived class CTORs
    int fwd_;

    // for sending messages in the forward direction (like a bcast)
    virtual int fwdsendto() = 0;
    virtual int fwdsend() = 0;
    virtual int fwdrecvfrom() = 0;
    virtual int fwdrecv() = 0;
  public:
    GlobalMsgIter(int nproc, int me, int root = 0);
    ~GlobalMsgIter();
    void backwards() { fwd_ = 0; i_ = n_ - 1; }
    void forwards() { fwd_ = 1; i_ = 0; }
    void next() { if (fwd_) i_++; else i_--; }
    int done() { return i_<0 || i_>=n_; }
    int n() { return n_; }
    int sendto() { return fwd_?fwdsendto():fwdrecvfrom(); }
    int send() { return fwd_?fwdsend():fwdrecv(); }
    int recvfrom() { return fwd_?fwdrecvfrom():fwdsendto(); }
    int recv() { return fwd_?fwdrecv():fwdsend(); }
};


class MessageGrp;
class MachineTopology: public DescribedClass {
  public:
    MachineTopology();
    MachineTopology(const Ref<KeyVal>&);
    ~MachineTopology();

    virtual Ref<GlobalMsgIter> global_msg_iter(const Ref<MessageGrp>&,
                                             int target) = 0;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
