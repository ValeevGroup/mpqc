//
// memrdma.h
// Based on memamsg.h
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

#ifndef _util_group_memrdma_h
#define _util_group_memrdma_h

#include <iostream>
#include <vector>

#include <util/group/pool.h>
#include <util/group/memmsg.h>

namespace sc {

/** The RDMAMemoryGrp abstract class specializes the MsgMemoryGrp
class.  It uses RDMA to implement global shared memory.  */
class RDMAMemoryGrp : public MsgMemoryGrp {
  protected:
    char *data_;

    virtual void retrieve_data(void *, int node, long offset, long size,
                               int lock) = 0;
    virtual void replace_data(void *, int node, long offset, long size,
                              int unlock) = 0;
    virtual void sum_data(double *data, int node, long doffset, long dsize) = 0;

    std::vector<Pool*> pools_;
    size_t default_pool_size_;
    void* malloc_region(size_t nbyte);
    void free_region(void*);
  public:
    RDMAMemoryGrp(const Ref<MessageGrp>& msg);
    RDMAMemoryGrp(const Ref<KeyVal>&);
    ~RDMAMemoryGrp();

    void *localdata();

    void set_localsize(size_t localsize);

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
