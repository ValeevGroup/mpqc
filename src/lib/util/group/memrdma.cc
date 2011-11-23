//
// memrdma.cc
// Based on memamsg.cc
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

#ifndef _util_group_memrdma_cc
#define _util_group_memrdma_cc

#include <util/class/scexception.h>
#include <util/misc/formio.h>
#include <util/group/pool.h>
#include <util/group/memrdma.h>
#include <util/group/memiter.h>
#include <stdexcept>

using namespace std;
using namespace sc;

#ifdef HAVE_HRECV
#  define DISABLE do { masktrap(1); ExEnv::outn().flush(); } while(0)
#  define ENABLE do { ExEnv::outn().flush(); masktrap(0); } while(0)
   extern "C" {
       long masktrap(long state);
     }
#else
#  define DISABLE ExEnv::outn().flush()
#  define ENABLE ExEnv::outn().flush()
#endif

#define PRINTF(args) do { DISABLE; \
                          ExEnv::outn() << scprintf args ; \
                          ExEnv::outn().flush(); \
                          ENABLE; \
                         } while(0)

#undef PRINTF
#define PRINTF(args)

///////////////////////////////////////////////////////////////////////
// Members for RDMAMemoryGrp

static ClassDesc RDMAMemoryGrp_cd(
  typeid(RDMAMemoryGrp),"RDMAMemoryGrp",1,"public MsgMemoryGrp",
  0, 0, 0);

RDMAMemoryGrp::RDMAMemoryGrp(const Ref<MessageGrp>& msg):
  MsgMemoryGrp(msg)
{
  data_ = 0;
  default_pool_size_ = 1000000;
}

RDMAMemoryGrp::RDMAMemoryGrp(const Ref<KeyVal>& keyval):
  MsgMemoryGrp(keyval)
{
  data_ = 0;
  default_pool_size_ = 1000000;
}

void*
RDMAMemoryGrp::malloc_region(size_t nbyte)
{
  void *data = 0;
  for (int i=0; data==0 && i<pools_.size(); i++) {
      data = pools_[i]->allocate(nbyte);
    }
  if (data == 0) {
      if (default_pool_size_ < nbyte) default_pool_size_ = nbyte * 2;
      else if (pools_.size() > 4) default_pool_size_ *= 2;
      void *pooldata = malloc_local(default_pool_size_);
      Pool *pool = new(pooldata) Pool(default_pool_size_);
      pools_.push_back(pool);
      data = pool->allocate(nbyte);
    }
  return data;
}

void
RDMAMemoryGrp::free_region(void*data)
{
  char *cdata = reinterpret_cast<char*>(data);
  for (int i=0; i<pools_.size(); i++) {
      char *pstart = reinterpret_cast<char*>(pools_[i]);
      if (cdata > pstart && cdata < &pstart[pools_[i]->size()]) {
          pools_[i]->release(data);
          return;
        }
    }
  throw ProgrammingError("could not find data to release in a Pool",
                         __FILE__, __LINE__, this->class_desc());
}

void
RDMAMemoryGrp::set_localsize(size_t localsize)
{
  for (int i=0; i<pools_.size(); i++) {
    void* pool_i_ptr = pools_[i];
    free_local(pool_i_ptr);
  }
  pools_.resize(0);

  MsgMemoryGrp::set_localsize(localsize);
}

void *
RDMAMemoryGrp::localdata()
{
  return data_;
}

RDMAMemoryGrp::~RDMAMemoryGrp()
{
  deactivate();
  //delete[] data_;
}

void *
RDMAMemoryGrp::obtain_writeonly(distsize_t offset, int size)
{
  void *data = malloc_region(size);
  return data;
}

void *
RDMAMemoryGrp::obtain_readwrite(distsize_t offset, int size)
{
  PRINTF(("RDMAMemoryGrp::obtain_readwrite entered\n"));
  void *data = malloc_region(size);
  MemoryIter i(data, offsets_, n());
  for (i.begin(offset, size); i.ready(); i.next()) {
      PRINTF(("RDMAMemoryGrp::obtain_readwrite: node = %d, "
              "int offset = %d, int size = %d\n",
              i.node(), i.offset()/sizeof(int), i.size()/sizeof(int)));
      retrieve_data(i.data(), i.node(), i.offset(), i.size(), 1);
    }
  PRINTF(("RDMAMemoryGrp::obtain_readwrite exiting\n"));
  return data;
}

void *
RDMAMemoryGrp::obtain_readonly(distsize_t offset, int size)
{
  void *data = malloc_region(size);
  PRINTF(("%d: RDMAMemoryGrp::obtain_readonly:"
          "overall: offset = %d size = %d\n",
          me(), offset, size));
  MemoryIter i(data, offsets_, n());
  for (i.begin(offset, size); i.ready(); i.next()) {
      PRINTF(("%d: RDMAMemoryGrp::obtain_readonly:working on:"
              "node = %d offset = %d size = %d\n",
              me(), i.node(), i.offset(), i.size()));
      PRINTF(("RDMAMemoryGrp::obtain_readonly: node = %d, "
              "int offset = %d, int size = %d\n",
              i.node(), i.offset()/sizeof(int), i.size()/sizeof(int)));
      retrieve_data(i.data(), i.node(), i.offset(), i.size(), 0);
    }
  return data;
}

void
RDMAMemoryGrp::sum_reduction(double *data, distsize_t doffset, int dsize)
{
  distsize_t offset = doffset * sizeof(double);
  int size = dsize * sizeof(double);
  MemoryIter i(data, offsets_, n());
  for (i.begin(offset, size); i.ready(); i.next()) {
      sum_data((double*)i.data(), i.node(), i.offset(), i.size());
    }
}

void
RDMAMemoryGrp::sum_reduction_on_node(double *data, size_t doffset,
                                          int dlength, int node)
{
  if (node == -1) node = me();

  sum_data(data, node, sizeof(double)*doffset, sizeof(double)*dlength);
}

void
RDMAMemoryGrp::release_readonly(void *data, distsize_t offset, int size)
{
  free_region(data);
}

void
RDMAMemoryGrp::release_writeonly(void *data, distsize_t offset, int size)
{
  MemoryIter i(data, offsets_, n());
  for (i.begin(offset, size); i.ready(); i.next()) {
      PRINTF(("RDMAMemoryGrp::release_write: node = %d, "
              "int offset = %d, int size = %d\n",
              i.node(), i.offset()/sizeof(int), i.size()/sizeof(int)));
      replace_data(i.data(), i.node(), i.offset(), i.size(), 0);
    }
  free_region(data);
}

void
RDMAMemoryGrp::release_readwrite(void *data, distsize_t offset, int size)
{
  MemoryIter i(data, offsets_, n());
  for (i.begin(offset, size); i.ready(); i.next()) {
      replace_data(i.data(), i.node(), i.offset(), i.size(), 1);
    }
  free_region(data);
}

void
RDMAMemoryGrp::print(ostream &o) const
{
  MemoryGrp::print(o);
}

#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
