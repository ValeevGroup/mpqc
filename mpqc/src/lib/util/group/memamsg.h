//
// memamsg.h
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

#ifndef _util_group_memamsg_h
#define _util_group_memamsg_h

#include <iostream>

#include <util/group/memmsg.h>

namespace sc {
    
class MemoryDataRequest {
  public:
    enum { NData = 6 };
    enum Request { Deactivate, Sync, Retrieve, Replace, DoubleSum };
  private:
    int data_[NData];
  public:
    MemoryDataRequest() {}
    MemoryDataRequest(Request r, int node = 0, int offset = 0, int size = 0,
                      int lock = 0, int serial = 0);
    void assign(Request r, int node, int offset, int size,
                int lock, int serial);
    void *data() const { return (void *) data_; }
    int nbytes() const { return sizeof(int)*NData; }

    const char *request_string() const;

    MemoryDataRequest::Request request() const { return (Request) data_[0]; }
    int node() const { return data_[1]; }
    int offset() const { return data_[2]; }
    int size() const { return data_[3]; }
    int serial_number() const { return data_[4]; }
    int lock() const { return data_[5]; }

    int touches_data() const {return request()!=Deactivate&&request()!=Sync;}

    // Sync messages only define one datum besides type and node
    int reactivate() const { return data_[2]; }

    void operator =(const MemoryDataRequest &r);

    void print(const char* msg = 0, std::ostream & o = ExEnv::out0()) const;
};

class MemoryDataRequestQueue {
  public:
    enum { MaxDepth = 1024 };
  private:
    MemoryDataRequest q_[MaxDepth];
    int n_;
  public:
    MemoryDataRequestQueue(): n_(0) {}
    int n() const { return n_; }
    void push(MemoryDataRequest&);
    void pop(MemoryDataRequest&);

    MemoryDataRequest& operator[](int i) { return q_[i]; }
    void clear() { n_ = 0; }
};

/** The ActiveMsgMemoryGrp abstract class specializes the MsgMemoryGrp
class.  It uses active messages to implement global shared memory.  */
class ActiveMsgMemoryGrp : public MsgMemoryGrp {
  protected:
    char *data_;

    virtual void retrieve_data(void *, int node, int offset, int size,
                               int lock) = 0;
    virtual void replace_data(void *, int node, int offset, int size,
                              int unlock) = 0;
    virtual void sum_data(double *data, int node, int doffset, int dsize) = 0;
  public:
    ActiveMsgMemoryGrp(const Ref<MessageGrp>& msg);
    ActiveMsgMemoryGrp(const Ref<KeyVal>&);
    ~ActiveMsgMemoryGrp();

    void set_localsize(size_t);
    void *localdata();

    void *obtain_writeonly(distsize_t offset, int size);
    void *obtain_readwrite(distsize_t offset, int size);
    void *obtain_readonly(distsize_t offset, int size);
    void release_readonly(void *data, distsize_t offset, int size);
    void release_writeonly(void *data, distsize_t offset, int size);
    void release_readwrite(void *data, distsize_t offset, int size);

    void sum_reduction(double *data, distsize_t doffset, int dsize);
    void sum_reduction_on_node(double *data, size_t doffset, int dsize,
                               int node = -1);

    void print(std::ostream &o = ExEnv::out0()) const;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
