//
// messpvm.h
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

#ifndef _util_group_messpvm_h
#define _util_group_messpvm_h

#include <util/group/message.h>

class PVMMessageGrp: public MessageGrp {
#define CLASSNAME PVMMessageGrp
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classd.h>
  protected:
    int *tids;

    int rtid;
    int rtag;
    int rlen;
  public:
    PVMMessageGrp();
    PVMMessageGrp(const RefKeyVal&);
    ~PVMMessageGrp();

    void raw_send(int target, void* data, int nbyte);
    void raw_recv(int sender, void* data, int nbyte);
    void raw_sendt(int target, int type, void* data, int nbyte);
    void raw_recvt(int type, void* data, int nbyte);

    int probet(int type);

    int last_source();
    int last_size();
    int last_type();
};

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
