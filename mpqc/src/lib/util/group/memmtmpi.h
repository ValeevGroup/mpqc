//
// memmtmpi.h
// based on memmpi.h
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

#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_memmtmpi_h
#define _util_group_memmtmpi_h

#include <fstream.h>

#include <util/group/message.h>
#include <util/group/memamsg.h>
#include <util/group/thread.h>

//. This MemoryGrp class requires a MT-safe MPI implementation.
//The default MessageGrp must be a MPIMessageGrp.  MPI must
//be safe with respect to the default ThreadGrp.  Alternately,
//a MessageGrp and a ThreadGrp can be passed to the constructor.
class MTMPIMemoryGrp: public ActiveMsgMemoryGrp {
#define CLASSNAME MTMPIMemoryGrp
#define HAVE_KEYVAL_CTOR
#include <util/class/classd.h>
  private:
    RefThreadGrp th_;

    int req_type_;
    int to_type_;
    int fr_type_;

    int active_;

    Thread *thread_;
    RefThreadLock mem_lock_;
    RefThreadLock print_lock_; // needed for debugging only
    ofstream hout; // handler out
    ofstream mout; // main thread out

    void init_mtmpimg();

    // parent class pure virtuals
    void retrieve_data(void *, int node, int offset, int size);
    void replace_data(void *, int node, int offset, int size);
    void sum_data(double *data, int node, int doffset, int dsize);

    friend class MTMPIThread;
  public:
    MTMPIMemoryGrp(const RefMessageGrp& msg, const RefThreadGrp &th);
    MTMPIMemoryGrp(const RefKeyVal &);
    ~MTMPIMemoryGrp();

    long lockcomm();
    void unlockcomm(long oldvalue);

    void activate();
    void deactivate();
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
