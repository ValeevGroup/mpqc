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
#include <util/group/pool.h>
#include <util/group/memamsg.h>
#include <util/group/memiter.h>

#ifdef HAVE_HRECV
#  define DISABLE do { masktrap(1); cout.flush(); } while(0)
#  define ENABLE do { cout.flush(); masktrap(0); } while(0)
   extern "C" {
       long masktrap(long state);
     }
#else
#  define DISABLE cout.flush()
#  define ENABLE cout.flush()
#endif

#define PRINTF(args) do { DISABLE; \
                          cout << scprintf args ; \
                          cout.flush(); \
                          ENABLE; \
                         } while(0)

#undef PRINTF
#define PRINTF(args)

///////////////////////////////////////////////////////////////////////
// The MemoryLockRequest class

MemoryLockRequest::MemoryLockRequest(Request r, int node,
                                     int start, int end)
{
  assign(r, node, start, end);
}

void
MemoryLockRequest::assign(Request r, int node,
                          int start, int end)
{
  data_[0] = (int) r;
  data_[1] = node;
  data_[2] = start;
  data_[3] = end;
}

///////////////////////////////////////////////////////////////////////
// The MemoryDataRequest class

static int request_serial_number;

MemoryDataRequest::MemoryDataRequest(Request r, int node,
                                     int offset, int size)
{
  assign(r, node, offset, size);
}

void
MemoryDataRequest::assign(Request r, int node,
                          int offset, int size)
{
  data_[0] = (int) r;
  data_[1] = node;
  data_[2] = offset;
  data_[3] = size;
  data_[4] = request_serial_number++;
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
      o << scprintf("%s \"%s\" offset = %5d, %5d bytes, %d-%d\n",
              msg, request_string(),
              offset(), size(), node(), serial_number());
    }
  o.flush();
}

void
MemoryDataRequest::operator =(const MemoryDataRequest &r)
{
  for (int i=0; i<NData; i++) {
      data_[i] = r.data_[i];
    }
}

///////////////////////////////////////////////////////////////////////
// The MemoryDataRequestQueue class

void
MemoryDataRequestQueue::push(MemoryDataRequest&r)
{
  if (n_ == MaxDepth) {
      cerr << scprintf("MemoryDataRequestQueue: MaxDepth exceeded\n");
      abort();
    }
  q_[n_] = r;
  n_++;
}

void
MemoryDataRequestQueue::pop(MemoryDataRequest&r)
{
  if (n_ == 0) {
      cerr << scprintf("MemoryDataRequestQueue: nothing to pop\n");
      abort();
    }
  n_--;
  r = q_[n_];
}

///////////////////////////////////////////////////////////////////////
// Members for ActiveMsgMemoryGrp

#define CLASSNAME ActiveMsgMemoryGrp
#define PARENTS public MsgMemoryGrp
#include <util/class/classia.h>
void *
ActiveMsgMemoryGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  MsgMemoryGrp::_castdown(cd);
  return do_castdowns(casts,cd);
}

ActiveMsgMemoryGrp::ActiveMsgMemoryGrp(const RefMessageGrp& msg):
  MsgMemoryGrp(msg)
{
  use_locks_for_reduction_ = 0;
  data_ = 0;
}

ActiveMsgMemoryGrp::ActiveMsgMemoryGrp(const RefKeyVal& keyval):
  MsgMemoryGrp(keyval)
{
  use_locks_for_reduction_ = 0;
  data_ = 0;
}

void
ActiveMsgMemoryGrp::set_localsize(int localsize)
{
  if (debug_) {
      cout << "ActiveMsgMemoryGrp::set_localsize(" << localsize << ")" << endl;
    }
  deactivate();
  MsgMemoryGrp::set_localsize(localsize);
  delete[] data_;
  data_ = new char[localsize];
  activate();
  if (debug_) {
      cout << "ActiveMsgMemoryGrp::set_localsize done: offsets:";
      for (int i=0; i<=n(); i++) {
          cout << " " << offset(i);
        }
      cout << endl;
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
  delete[] data_;
}

void *
ActiveMsgMemoryGrp::obtain_writeonly(distsize_t offset, int size)
{
  if (use_locks_) {
      send_lock_request(MemoryLockRequest::WriteOnly, offset, size);
      wait_for_lock();
    }
  void *data = (void *) new char[size];
  return data;
}

void *
ActiveMsgMemoryGrp::obtain_readwrite(distsize_t offset, int size)
{
  PRINTF(("ActiveMsgMemoryGrp::obtain_readwrite entered\n"));
  if (use_locks_) {
      send_lock_request(MemoryLockRequest::ReadWrite, offset, size);
      wait_for_lock();
    }
  void *data = (void *) new char[size];
  MemoryIter i(data, offsets_, n());
  for (i.begin(offset, size); i.ready(); i.next()) {
      if (i.node() == me()) {
          PRINTF(("ActiveMsgMemoryGrp::obtain_readwrite: local copy\n"));
          memcpy(i.data(), &data_[i.offset()], i.size());
        }
      else {
          PRINTF(("ActiveMsgMemoryGrp::obtain_readwrite: node = %d, "
                  "int offset = %d, int size = %d\n",
                  i.node(), i.offset()/sizeof(int), i.size()/sizeof(int)));
          retrieve_data(i.data(), i.node(), i.offset(), i.size());
        }
    }
  PRINTF(("ActiveMsgMemoryGrp::obtain_readwrite exiting\n"));
  return data;
}

void *
ActiveMsgMemoryGrp::obtain_readonly(distsize_t offset, int size)
{
  if (use_locks_) {
      send_lock_request(MemoryLockRequest::ReadOnly, offset, size);
      wait_for_lock();
    }
  void *data = (void *) new char[size];
  PRINTF(("%d: ActiveMsgMemoryGrp::obtain_readonly:"
          "overall: offset = %d size = %d\n",
          me(), offset, size));
  MemoryIter i(data, offsets_, n());
  for (i.begin(offset, size); i.ready(); i.next()) {
      PRINTF(("%d: ActiveMsgMemoryGrp::obtain_readonly:working on:"
              "node = %d offset = %d size = %d\n",
              me(), i.node(), i.offset(), i.size()));
      if (i.node() == me()) {
          PRINTF(("%d: ActiveMsgMemoryGrp::obtain_readonly: local: "
                  "offset = %d size = %d\n", me(), i.offset(), i.size()));
          memcpy(i.data(), &data_[i.offset()], i.size());
        }
      else {
          PRINTF(("ActiveMsgMemoryGrp::obtain_readonly: node = %d, "
                  "int offset = %d, int size = %d\n",
                  i.node(), i.offset()/sizeof(int), i.size()/sizeof(int)));
          retrieve_data(i.data(), i.node(), i.offset(), i.size());
        }
    }
  return data;
}

void
ActiveMsgMemoryGrp::sum_reduction(double *data, distsize_t doffset, int dsize)
{
  distsize_t offset = doffset * sizeof(double);
  int size = dsize * sizeof(double);
  // Locks are usually implicit, assuming that only one active message
  // handler is active at a time.
  if (use_locks_for_reduction_) {
      send_lock_request(MemoryLockRequest::Reduce, offset, size);
      wait_for_lock();
    }
  MemoryIter i(data, offsets_, n());
  for (i.begin(offset, size); i.ready(); i.next()) {
      if (i.node() == me()) {
          int chunkdsize = i.size()/sizeof(double);
          double *chunkdata = (double*) &data_[i.offset()];
          double *tmp = (double*) i.data();
          PRINTF(("%d: summing %d doubles from 0x%x to 0x%x\n",
                  me(), chunkdsize, tmp, chunkdata));
          long oldlock = lockcomm();
          for (int j=0; j<chunkdsize; j++) {
              *chunkdata++ += *tmp++;
            }
          unlockcomm(oldlock);
        }
      else {
          sum_data((double*)i.data(), i.node(), i.offset(), i.size());
        }
    }
  if (use_locks_for_reduction_) {
      send_lock_request(MemoryLockRequest::RelReduce, offset, size);
    }
}

void
ActiveMsgMemoryGrp::sum_reduction_on_node(double *data, int doffset,
                                          int dlength, int node)
{
  if (node == -1) node = me();

  if (node == me()) {
      double *localdata = (double*) &data_[sizeof(double)*doffset];
      long oldlock = lockcomm();
      for (int j=0; j<dlength; j++) {
          *localdata++ += *data++;
        }
      unlockcomm(oldlock);
    }
  else {
      sum_data(data, node, sizeof(double)*doffset, sizeof(double)*dlength);
    }
}

void
ActiveMsgMemoryGrp::release_read(void *data, distsize_t offset, int size)
{
  if (use_locks_) {
      send_lock_request(MemoryLockRequest::RelRead, offset, size);
    }
  delete[] (char*) data;
}

void
ActiveMsgMemoryGrp::release_write(void *data, distsize_t offset, int size)
{
  if (use_locks_) {
      send_lock_request(MemoryLockRequest::RelWrite, offset, size);
    }
  MemoryIter i(data, offsets_, n());
  for (i.begin(offset, size); i.ready(); i.next()) {
      if (i.node() == me()) {
          PRINTF(("ActiveMsgMemoryGrp::release_write: local\n"));
          PRINTF(("  i.offset() = %d i.data() = 0x%x i.size() = %d\n",
                  i.offset(), i.data(), i.size()));
          PRINTF(("  &data_[i.offset()] = 0x%x\n", &data_[i.offset()]));
          memcpy(&data_[i.offset()], i.data(), i.size());
        }
      else {
          PRINTF(("ActiveMsgMemoryGrp::release_write: node = %d, "
                  "int offset = %d, int size = %d\n",
                  i.node(), i.offset()/sizeof(int), i.size()/sizeof(int)));
          replace_data(i.data(), i.node(), i.offset(), i.size());
        }
    }
  delete[] (char*) data;
}

void
ActiveMsgMemoryGrp::print(ostream &o) const
{
  MemoryGrp::print(o);
}

void
ActiveMsgMemoryGrp::send_lock_request(MemoryLockRequest::Request req,
                                      distsize_t offset, int size)
{
  cerr << scprintf("%d: %s: cannot use memory locks\n", me(), class_name());
  abort();
}

void
ActiveMsgMemoryGrp::wait_for_lock()
{
  cerr << scprintf("%d: %s: cannot use memory locks\n", me(), class_name());
  abort();
}

long
ActiveMsgMemoryGrp::lockcomm()
{
    return 0;
}

void
ActiveMsgMemoryGrp::unlockcomm(long oldvalue)
{
}

#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
