//
// memmid.h
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

#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_memmid_h
#define _util_group_memmid_h

#include <util/group/memamsg.h>

// This is used for memory handler that use message identifiers to
// keep trace of messages.

class MIDMemoryGrp: public ActiveMsgMemoryGrp {
#define CLASSNAME MIDMemoryGrp
#include <util/class/classd.h>
  public:
    // This is public so memory handler functions can call it.
    void handler(long *mid = 0);
  protected:
    void handler(MemoryDataRequest&, long *mid = 0);
    
    void retrieve_data(void *, int node, int offset, int size);
    void replace_data(void *, int node, int offset, int size);
    void sum_data(double *data, int node, int doffset, int dsize);

    int active_;

    int data_request_type_;
    int data_type_to_handler_;
    int data_type_from_handler_;

    long data_request_mid_;

    MemoryDataRequest data_request_buffer_;

    int nsync_;

    int use_acknowledgments_;
    int use_active_messages_;

    void print_memreq(MemoryDataRequest &req,
                      const char * = 0, int target = -1);

    void do_wait(const char *msg, int mid,
                 MemoryDataRequestQueue &q, size_t expectedsize,
                 int node = -1 /*not needed except for debugging*/);
    void flush_queue(MemoryDataRequestQueue &q);

    virtual long lockcomm() = 0;
    virtual void unlockcomm(long oldvalue) = 0;
    virtual long send(void* data, int nbytes, int node, int type) = 0;
    virtual long recv(void* data, int nbytes, int node, int type) = 0;
    virtual long postrecv(void *data, int nbytes, int type) = 0;
    virtual long wait(long, long = -1) = 0;
    virtual int probe(long);

    virtual void got_data_request_mid();
  public:
    MemoryDataRequest &data_request_buffer() { return data_request_buffer_; }

    MIDMemoryGrp(const RefMessageGrp& msg);
    MIDMemoryGrp(const RefKeyVal& msg);
    ~MIDMemoryGrp();

    void activate();
    void deactivate();

    void sync();
    void catchup();
};

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
