//
// memmtmpi.cc
// based on memmpi.cc
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

#ifndef _util_group_memmtmpi_cc
#define _util_group_memmtmpi_cc

#ifdef __GNUC__
#pragma implementation
#endif

#include <assert.h>

#include <util/misc/formio.h>
#include <util/group/messmpi.h>
#include <util/group/memmtmpi.h>

#include <mpi.h>

using namespace std;

static const int dbufsize = 32768;

///////////////////////////////////////////////////////////////////////
// The MTMPIThread class

class MTMPIThread: public Thread {
  private:
    MTMPIMemoryGrp *mem_;
    int req_tag_;
    int tag_;
    double chunk[dbufsize];
  public:
    MTMPIThread(MTMPIMemoryGrp *, int reqtype, int tag);
    void run();
};

MTMPIThread::MTMPIThread(MTMPIMemoryGrp *mem,
                         int reqtype, int tag)
{
  mem_ = mem;
  req_tag_ = reqtype;
  tag_ = tag;
}

void
MTMPIThread::run()
{
  int i;
  int dsize;
  int dremain;
  int doffset;
  long l;
  MemoryDataRequest req;
  MPI_Status status;
  while (1) {
      MPI_Recv(req.data(),req.nbytes(),MPI_BYTE,MPI_ANY_SOURCE,
               req_tag_,mem_->comm_,&status);
      int rtag = req.serial_number();
      if (mem_->debug_) {
          mem_->print_lock_->lock();
          req.print("RECV",mem_->hout);
          mem_->print_lock_->unlock();
        }
      if (req.touches_data()) {
          assert(req.size() >= 0);
          if (req.offset()+req.size() > mem_->localsize()) {
              mem_->print_lock_->lock();
              req.print("BAD RECV");
              ExEnv::out() << "mem_->localsize() = " << mem_->localsize() << endl;
              mem_->print_lock_->lock();
            }
          assert(req.offset()+req.size() <= mem_->localsize());
        }
      switch (req.request()) {
      case MemoryDataRequest::Deactivate:
          return;
      case MemoryDataRequest::Retrieve:
          if (req.lock())
              mem_->obtain_local_lock(req.offset(), req.offset()+req.size());
          MPI_Send(&mem_->data_[req.offset()],req.size(),MPI_BYTE,
                   req.node(),rtag,mem_->comm_);
          break;
      case MemoryDataRequest::Replace:
          MPI_Send(&tag_,1,MPI_INTEGER,req.node(),rtag,mem_->comm_);
          MPI_Recv(&mem_->data_[req.offset()],req.size(),MPI_BYTE,
                   req.node(),tag_,mem_->comm_,&status);
          if (req.lock())
              mem_->release_local_lock(req.offset(), req.offset()+req.size());
          break;
      case MemoryDataRequest::DoubleSum:
          MPI_Send(&tag_,1,MPI_INTEGER,req.node(),rtag,mem_->comm_);
          dsize = req.size()/sizeof(double);
          dremain = dsize;
          doffset = req.offset()/sizeof(double);
          mem_->obtain_local_lock(req.offset(), req.offset()+req.size());
          while(dremain>0) {
              int dchunksize = dbufsize;
              if (dremain < dchunksize) dchunksize = dremain;
              MPI_Recv(chunk,dchunksize,MPI_DOUBLE,
                       req.node(),tag_,mem_->comm_,&status);
              double *source_data = &((double*)mem_->data_)[doffset];
              for (i=0; i<dchunksize; i++) {
                  source_data[i] += chunk[i];
                }
              dremain -= dchunksize;
              doffset += dchunksize;
            }
          mem_->release_local_lock(req.offset(), req.offset()+req.size());
          break;
      default:
          mem_->print_lock_->lock();
          ExEnv::out() << "MTMPIThread: bad memory data request" << endl;
          mem_->print_lock_->unlock();
          abort();
        }
    }
}

///////////////////////////////////////////////////////////////////////
// The MTMPIMemoryGrp class

static ClassDesc MTMPIMemoryGrp_cd(
  typeid(MTMPIMemoryGrp),"MTMPIMemoryGrp",1,"public ActiveMsgMemoryGrp",
  0, create<MTMPIMemoryGrp>, 0);

MTMPIMemoryGrp::MTMPIMemoryGrp(const Ref<MessageGrp>& msg,
                               const Ref<ThreadGrp>& th):
  ActiveMsgMemoryGrp(msg)
{
  if (debug_) ExEnv::out() << "MTMPIMemoryGrp CTOR entered" << endl;

  th_ = th;

  init_mtmpimg(th_->nthread());
}

MTMPIMemoryGrp::MTMPIMemoryGrp(const Ref<KeyVal>& keyval):
  ActiveMsgMemoryGrp(keyval)
{
  if (debug_) ExEnv::out() << "MTMPIMemoryGrp keyval CTOR entered" << endl;

  th_ = ThreadGrp::get_default_threadgrp();

  KeyValValueint nthreaddef(th_->nthread());
  int nthread = keyval->intvalue("num_threads",nthreaddef);
  ExEnv::out() << node0 << indent << "MTMPIMemoryGrp: num_threads = " << nthread << endl;

  init_mtmpimg(nthread);
}

MTMPIMemoryGrp::~MTMPIMemoryGrp()
{
  deactivate();
  for (int i=0; i<th_->nthread()-1; i++) {
      delete thread_[i];
    }
  delete[] thread_;
}

void
MTMPIMemoryGrp::init_mtmpimg(int nthread)
{
  int i;
  active_ = 0;

  if (nthread < 2) nthread = 2;
  th_ = th_->clone(nthread);
  nthread = th_->nthread();
  if (nthread < 2) {
      ExEnv::out() << "MTMPIMemoryGrp didn't get enough threads" << endl;
      abort();
    }

  if (debug_) {
      char name[256];
      sprintf(name, "mpqc.hand.%d", me());
      hout.open(name);
      sprintf(name, "mpqc.main.%d", me());
      mout.open(name);
    }

  MPI_Comm_dup(MPI_COMM_WORLD, &comm_);

  serial_ = 0;
  req_tag_ = 15001;

  serial_lock_ = th_->new_lock();

  thread_ = new Thread*[nthread-1];
  th_->add_thread(0,0);
  for (i=1; i<nthread; i++) {
      thread_[i-1] = new MTMPIThread(this,req_tag_,req_tag_ + i);
      th_->add_thread(i,thread_[i-1]);
    }
  print_lock_ = th_->new_lock();
}

int
MTMPIMemoryGrp::serial()
{
  serial_lock_->lock();
  int r = serial_;
  serial_++;
  if (serial_ == req_tag_) serial_ = 0;
  serial_lock_->unlock();
  return r;
}

void
MTMPIMemoryGrp::retrieve_data(void *data, int node, int offset, int size,
                              int lock)
{
  MemoryDataRequest req(MemoryDataRequest::Retrieve,me(),offset,size,lock,
                        serial());
  int tag = req.serial_number();

  // send the request
  if (debug_) {
      print_lock_->lock();
      req.print("SEND",mout);
      print_lock_->unlock();
    }
  MPI_Send(req.data(),req.nbytes(),MPI_BYTE,node,req_tag_,comm_);

  // receive the data
  MPI_Status status;
  MPI_Recv(data,size,MPI_BYTE,node,tag,comm_,&status);
}

void
MTMPIMemoryGrp::replace_data(void *data, int node, int offset, int size,
                             int unlock)
{
  MemoryDataRequest req(MemoryDataRequest::Replace,me(),offset,size,unlock,
                        serial());
  int tag = req.serial_number();

  if (debug_) {
      print_lock_->lock();
      req.print("SEND",mout);
      print_lock_->unlock();
    }
  MPI_Send(req.data(),req.nbytes(),MPI_BYTE,node,req_tag_,comm_);

  // wait for the go ahead message
  MPI_Status status;
  int rtag;
  MPI_Recv(&rtag,1,MPI_INTEGER,node,tag,comm_,&status);

  // send the data
  MPI_Send(data,size,MPI_BYTE,node,rtag,comm_);
}

void
MTMPIMemoryGrp::sum_data(double *data, int node, int offset, int size)
{
  MemoryDataRequest req(MemoryDataRequest::DoubleSum,me(),offset,size,
                        0, serial());
  int tag = req.serial_number();

  // send the request
  if (debug_) {
      print_lock_->lock();
      req.print("SEND",mout);
      print_lock_->unlock();
    }
  MPI_Send(req.data(),req.nbytes(),MPI_BYTE,node,req_tag_,comm_);

  // wait for the go ahead message
  MPI_Status status;
  int rtag;
  MPI_Recv(&rtag,1,MPI_INTEGER,node,tag,comm_,&status);

  int dsize = size/sizeof(double);
  int dremain = dsize;
  int dcurrent = 0;
  while(dremain>0) {
      int dchunksize = dbufsize;
      if (dremain < dchunksize) dchunksize = dremain;
      // send the data
      MPI_Send(&data[dcurrent],dchunksize,MPI_DOUBLE,
               node,rtag,comm_);
      dcurrent += dchunksize;
      dremain -= dchunksize;
    }
}

void
MTMPIMemoryGrp::activate()
{
  if (active_) return;
  active_ = 1;

  th_->start_threads();
}

void
MTMPIMemoryGrp::deactivate()
{
  if (!active_) return;
  active_ = 0;

  // send a shutdown message
  MemoryDataRequest req(MemoryDataRequest::Deactivate);
  if (debug_) {
      print_lock_->lock();
      req.print("SEND",mout);
      print_lock_->unlock();
    }
  for (int i=1; i<th_->nthread(); i++) {
      MPI_Send(req.data(),req.nbytes(),MPI_BYTE,me(),req_tag_,comm_);
    }

  // wait on the thread to shutdown
  th_->wait_threads();
}

void
MTMPIMemoryGrp::sync()
{
  // This method for sync eliminates the need for an "all done"
  // message after a memory operation completes.  However, a
  // handshake at the beginning of memory transactions are needed
  // so msg_->sync() complete before all remote memory operations
  // have at least begun.  Otherwise, I could deactivate before
  // all messages arrive.
  msg_->sync();
  if (active_) {
      deactivate();
      // At this point all memory operations are complete;
      // start up for the next phase.
      activate();
    }
}

#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
