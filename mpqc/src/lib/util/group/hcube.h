//
// hcube.h
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

#ifndef _util_group_hcube_h
#define _util_group_hcube_h

#include <util/group/topology.h>

namespace sc {

class HypercubeGMI: public GlobalMsgIter {
  private:
    int nhyper_;
  protected:
    int fwdsendto();
    int fwdsend();
    int fwdrecvfrom();
    int fwdrecv();
  public:
    HypercubeGMI(int nproc, int me, int root);
    ~HypercubeGMI();
};

// This utilitizes a hypercube topology, but will work for any number of
// nodes.
class HypercubeTopology: public MachineTopology {
  public:
    HypercubeTopology();
    HypercubeTopology(const Ref<KeyVal>&);
    ~HypercubeTopology();
    Ref<GlobalMsgIter> global_msg_iter(const Ref<MessageGrp>&, int target);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
