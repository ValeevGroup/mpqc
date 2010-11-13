//
// memamsg.cc
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

#ifndef _util_group_memamsg_cc
#define _util_group_memamsg_cc

#ifdef __GNUC__
#pragma implementation
#endif

#include <util/misc/formio.h>
#include <util/misc/consumableresources.h>
#include <util/group/pool.h>
#include <util/group/memamsg.h>
#include <util/group/memiter.h>

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

// comment out these two lines to produce diagnostic output
#undef PRINTF
#define PRINTF(args)

///////////////////////////////////////////////////////////////////////
// The MemoryDataRequest class

MemoryDataRequest::MemoryDataRequest(Request r, int node,
                                     long offset, long size, int lock,
                                     int serial_number)
{
  assign(r, node, offset, size, lock, serial_number);
}

void
MemoryDataRequest::assign(Request r, int node,
                          long offset, long size, int lock,
                          int serial_number)
{
  data_.request = r;
  data_.node = node;
  data_.offset = offset;
  data_.size = size;
  data_.serial_number = serial_number;
  data_.lock = lock;
}

const char *
MemoryDataRequest::request_string() const
{
  switch (request()) {
  case MemoryDataRequest::Deactivate:
      return "Deactivate";
  case MemoryDataRequest::Retrieve:
      return "Retrieve";
  case MemoryDataRequest::Replace:
      return "Replace";
  case MemoryDataRequest::DoubleSum:
      return "DoubleSum";
  case MemoryDataRequest::Sync:
      return "Sync";
  default:
      return "BadRequest";
    }
}

void
MemoryDataRequest::print(const char *msg, ostream & o) const
{
  if (msg == 0) msg = "";

  o.flush();
  if (request() == Sync) {
      o << scprintf("%s \"%s\" %d-%d reactivate = %d\n",
              msg, request_string(), node(), serial_number(), reactivate());
    }
  else {
      o << scprintf("%s \"%s\" offset = %5d, %5d bytes, %d-%d, %s\n",
              msg, request_string(),
              offset(), size(), node(), serial_number(),
              (lock()?"lock":"nolock"));
    }
  o.flush();
}

void
MemoryDataRequest::operator =(const MemoryDataRequest &r)
{
  data_ = r.data_;
}

///////////////////////////////////////////////////////////////////////
// The MemoryDataRequestQueue class

void
MemoryDataRequestQueue::push(MemoryDataRequest&r)
{
  if (n_ == MaxDepth) {
      ExEnv::errn() << scprintf("MemoryDataRequestQueue: MaxDepth exceeded\n");
      abort();
    }
  q_[n_] = r;
  n_++;
}

void
MemoryDataRequestQueue::pop(MemoryDataRequest&r)
{
  if (n_ == 0) {
      ExEnv::errn() << scprintf("MemoryDataRequestQueue: nothing to pop\n");
      abort();
    }
  n_--;
  r = q_[n_];
}

///////////////////////////////////////////////////////////////////////
// Members for ActiveMsgMemoryGrp

static ClassDesc ActiveMsgMemoryGrp_cd(
  typeid(ActiveMsgMemoryGrp),"ActiveMsgMemoryGrp",1,"public MsgMemoryGrp",
  0, 0, 0);

ActiveMsgMemoryGrp::ActiveMsgMemoryGrp(const Ref<MessageGrp>& msg):
  MsgMemoryGrp(msg)
{
  data_ = 0;
}

ActiveMsgMemoryGrp::ActiveMsgMemoryGrp(const Ref<KeyVal>& keyval):
  MsgMemoryGrp(keyval)
{
  data_ = 0;
}

void
ActiveMsgMemoryGrp::set_localsize(size_t localsize)
{
  if (debug_) {
      ExEnv::out0() << "ActiveMsgMemoryGrp::set_localsize(" << localsize << ")" << endl;
    }
  deactivate();
  MsgMemoryGrp::set_localsize(localsize);
  deallocate(data_);
  data_ = allocate<char>(localsize);
  activate();
  if (debug_) {
      ExEnv::out0() << "ActiveMsgMemoryGrp::set_localsize done: offsets:";
      for (int i=0; i<=n(); i++) {
          ExEnv::out0() << " " << double(offset(i));
        }
      ExEnv::out0() << endl;
    }
}

void *
ActiveMsgMemoryGrp::localdata()
{
  return data_;
}

ActiveMsgMemoryGrp::~ActiveMsgMemoryGrp()
{
  deactivate();
  deallocate(data_);
}

void *
ActiveMsgMemoryGrp::obtain_writeonly(distsize_t offset, int size)
{
  void *data = this->malloc_local(size);
  return data;
}

void *
ActiveMsgMemoryGrp::obtain_readwrite(distsize_t offset, int size)
{
  PRINTF(("%d: entering ActiveMsgMemoryGrp::obtain_readwrite: overall: offset = %ld size = %d\n",
          me(), (size_t)offset, size));
  void *data = this->malloc_local(size);
  MemoryIter i(data, offsets_, n());
  for (i.begin(offset, size); i.ready(); i.next()) {
    PRINTF(("%d: ActiveMsgMemoryGrp::obtain_readwrite: working on:"
            "node = %d offset = %d size = %d\n",
            me(), i.node(), i.offset(), i.size()));
      if (i.node() == me()) {
          PRINTF(("%d: ActiveMsgMemoryGrp::obtain_readwrite: local copy: "
                  "offset = %d size = %d\n", me(), i.offset(), i.size()));
          obtain_local_lock(i.offset(), i.offset()+i.size());
          memcpy(i.data(), &data_[i.offset()], i.size());
        }
      else {
          PRINTF(("%d: ActiveMsgMemoryGrp::obtain_readwrite: node = %d, "
                  "int offset = %d, int size = %d\n",
                  me(), i.node(), i.offset()/sizeof(int), i.size()/sizeof(int)));
          retrieve_data(i.data(), i.node(), i.offset(), i.size(), 1);
        }
    }
  PRINTF(("%d: exiting ActiveMsgMemoryGrp::obtain_readwrite\n", me()));
  return data;
}

void *
ActiveMsgMemoryGrp::obtain_readonly(distsize_t offset, int size)
{
  void *data = this->malloc_local(size);
  PRINTF(("%d: entering ActiveMsgMemoryGrp::obtain_readonly:"
          "overall: offset = %ld size = %d\n",
          me(), (size_t)offset, size));
  MemoryIter i(data, offsets_, n());
  for (i.begin(offset, size); i.ready(); i.next()) {
      PRINTF(("%d: ActiveMsgMemoryGrp::obtain_readonly: working on:"
              "node = %d offset = %d size = %d\n",
              me(), i.node(), i.offset(), i.size()));
      if (i.node() == me()) {
          PRINTF(("%d: ActiveMsgMemoryGrp::obtain_readonly: local: "
                  "offset = %d size = %d\n", me(), i.offset(), i.size()));
          memcpy(i.data(), &data_[i.offset()], i.size());
        }
      else {
          PRINTF(("%d: ActiveMsgMemoryGrp::obtain_readonly: node = %d, "
                  "int offset = %d, int size = %d\n",
                  me(), i.node(), i.offset()/sizeof(int), i.size()/sizeof(int)));
          retrieve_data(i.data(), i.node(), i.offset(), i.size(), 0);
        }
    }
  PRINTF(("%d: exiting ActiveMsgMemoryGrp::obtain_readonly\n", me()));
  return data;
}

void
ActiveMsgMemoryGrp::sum_reduction(double *data, distsize_t doffset, int dsize)
{
  distsize_t offset = doffset * sizeof(double);
  int size = dsize * sizeof(double);
  MemoryIter i(data, offsets_, n());
  for (i.begin(offset, size); i.ready(); i.next()) {
      if (i.node() == me()) {
          int chunkdsize = i.size()/sizeof(double);
          double *chunkdata = (double*) &data_[i.offset()];
          double *tmp = (double*) i.data();
          PRINTF(("%d: summing %d doubles from 0x%x to 0x%x\n",
                  me(), chunkdsize, tmp, chunkdata));
          obtain_local_lock(i.offset(), i.offset()+i.size());
          std::copy(tmp, tmp+chunkdsize, chunkdata);
#if 0
          for (int j=0; j<chunkdsize; j++) {
              *chunkdata++ += *tmp++;
            }
#endif
          release_local_lock(i.offset(), i.offset()+i.size());
        }
      else {
          sum_data((double*)i.data(), i.node(), i.offset(), i.size());
        }
    }
}

void
ActiveMsgMemoryGrp::sum_reduction_on_node(double *data, size_t doffset,
                                          int dlength, int node)
{
  if (node == -1) node = me();

  if (node == me()) {
      double *localdata = (double*) &data_[sizeof(double)*doffset];
      obtain_local_lock(sizeof(double)*doffset,
                        sizeof(double)*(doffset+dlength));
      for (int j=0; j<dlength; j++) {
          *localdata++ += *data++;
        }
      release_local_lock(sizeof(double)*doffset,
                         sizeof(double)*(doffset+dlength));
    }
  else {
      sum_data(data, node, sizeof(double)*doffset, sizeof(double)*dlength);
    }
}

void
ActiveMsgMemoryGrp::release_readonly(void *data, distsize_t offset, int size)
{
  this->free_local(data);
}

void
ActiveMsgMemoryGrp::release_writeonly(void *data, distsize_t offset, int size)
{
  MemoryIter i(data, offsets_, n());
  for (i.begin(offset, size); i.ready(); i.next()) {
      if (i.node() == me()) {
          PRINTF(("ActiveMsgMemoryGrp::release_write: local\n"));
          PRINTF(("  i.offset() = %d i.data() = 0x%x i.size() = %d\n",
                  i.offset(), i.data(), i.size()));
          PRINTF(("  &data_[i.offset()] = 0x%x\n", &data_[i.offset()]));
          const double* src_start = static_cast<const double*> (i.data());
          const double* src_fence = src_start + i.size();
          std::copy(src_start, src_fence, &data_[i.offset()]);
          //memcpy(&data_[i.offset()], i.data(), i.size());
        }
      else {
          PRINTF(("ActiveMsgMemoryGrp::release_write: node = %d, "
                  "int offset = %d, int size = %d\n",
                  i.node(), i.offset()/sizeof(int), i.size()/sizeof(int)));
          replace_data(i.data(), i.node(), i.offset(), i.size(), 0);
        }
    }
  this->free_local(data);
}

void
ActiveMsgMemoryGrp::release_readwrite(void *data, distsize_t offset, int size)
{
  PRINTF(("%d: entering ActiveMsgMemoryGrp::release_readwrite:"
          "overall: offset = %ld size = %d\n",
          me(), (size_t)offset, size));
  MemoryIter i(data, offsets_, n());
  for (i.begin(offset, size); i.ready(); i.next()) {
    PRINTF(("%d: ActiveMsgMemoryGrp::release_readwrite: working on:"
            "node = %d offset = %d size = %d\n",
            me(), i.node(), i.offset(), i.size()));
      if (i.node() == me()) {
        PRINTF(("%d: ActiveMsgMemoryGrp::release_readwrite: local copy: "
                "offset = %d size = %d\n", me(), i.offset(), i.size()));
        const double* src_start = static_cast<const double*> (i.data());
        const double* src_fence = src_start + i.size();
        std::copy(src_start, src_fence, &data_[i.offset()]);
          //memcpy(&data_[i.offset()], i.data(), i.size());
          release_local_lock(i.offset(), i.offset()+i.size());
        }
      else {
        PRINTF(("%d: ActiveMsgMemoryGrp::release_readwrite: node = %d, "
                "int offset = %d, int size = %d\n",
                me(), i.node(), i.offset()/sizeof(int), i.size()/sizeof(int)));
          replace_data(i.data(), i.node(), i.offset(), i.size(), 1);
        }
    }
  this->free_local(data);
  PRINTF(("%d: exiting ActiveMsgMemoryGrp::release_readwrite\n", me()));
}

void
ActiveMsgMemoryGrp::print(ostream &o) const
{
  MemoryGrp::print(o);
}

#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
