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

static const int dbufsize = 32768;

///////////////////////////////////////////////////////////////////////
// The MTMPIThread class

class MTMPIThread: public Thread {
  private:
    MTMPIMemoryGrp *mem_;
    int req_type_;
    int to_type_;
    int fr_type_;
  public:
    MTMPIThread(MTMPIMemoryGrp *, int reqtype, int totype, int fromtype);
    void run();
};

MTMPIThread::MTMPIThread(MTMPIMemoryGrp *mem,
                         int reqtype, int totype, int fromtype)
{
  mem_ = mem;
  req_type_ = reqtype;
  to_type_ = totype;
  fr_type_ = fromtype;
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
  double chunk[dbufsize];
  char junk;
  while (1) {
      MPI_Recv(req.data(),req.nbytes(),MPI_BYTE,MPI_ANY_SOURCE,
               req_type_,MPI_COMM_WORLD,&status);
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
              cout << "mem_->localsize() = " << mem_->localsize() << endl;
              mem_->print_lock_->lock();
            }
          assert(req.offset()+req.size() <= mem_->localsize());
        }
      switch (req.request()) {
      case MemoryDataRequest::Deactivate:
          return;
      case MemoryDataRequest::Retrieve:
          MPI_Send(&mem_->data_[req.offset()],req.size(),MPI_BYTE,
                   req.node(),fr_type_,MPI_COMM_WORLD);
          break;
      case MemoryDataRequest::Replace:
          MPI_Send(&junk,1,MPI_BYTE,req.node(),fr_type_,MPI_COMM_WORLD);
          MPI_Recv(&mem_->data_[req.offset()],req.size(),MPI_BYTE,
                   req.node(),to_type_,MPI_COMM_WORLD,&status);
          MPI_Send(&junk,1,MPI_BYTE,req.node(),fr_type_,MPI_COMM_WORLD);
          break;
      case MemoryDataRequest::DoubleSum:
          MPI_Send(&junk,1,MPI_BYTE,req.node(),fr_type_,MPI_COMM_WORLD);
          dsize = req.size()/sizeof(double);
          dremain = dsize;
          doffset = req.offset()/sizeof(double);
          l = mem_->lockcomm();
          while(dremain>0) {
              int dchunksize = dbufsize;
              if (dremain < dchunksize) dchunksize = dremain;
              MPI_Recv(chunk,dchunksize,MPI_DOUBLE,
                       req.node(),to_type_,MPI_COMM_WORLD,&status);
              double *source_data = &((double*)mem_->data_)[doffset];
              for (i=0; i<dchunksize; i++) {
                  source_data[i] += chunk[i];
                }
              dremain -= dchunksize;
              doffset += dchunksize;
            }
          mem_->unlockcomm(l);
          MPI_Send(&junk,1,MPI_BYTE,req.node(),fr_type_,MPI_COMM_WORLD);
          break;
      default:
          mem_->print_lock_->lock();
          cout << "MTMPIThread: bad memory data request" << endl;
          mem_->print_lock_->unlock();
          abort();
        }
    }
}

///////////////////////////////////////////////////////////////////////
// The MTMPIMemoryGrp class

#define CLASSNAME MTMPIMemoryGrp
#define HAVE_KEYVAL_CTOR
#define PARENTS public ActiveMsgMemoryGrp
#include <util/class/classi.h>
void *
MTMPIMemoryGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  ActiveMsgMemoryGrp::_castdown(cd);
  return do_castdowns(casts,cd);
}

MTMPIMemoryGrp::MTMPIMemoryGrp(const RefMessageGrp& msg,
                               const RefThreadGrp& th):
  ActiveMsgMemoryGrp(msg)
{
  if (debug_) cout << "MTMPIMemoryGrp CTOR entered" << endl;

  th_ = th;

  init_mtmpimg(th_->nthread() + 1);
}

MTMPIMemoryGrp::MTMPIMemoryGrp(const RefKeyVal& keyval):
  ActiveMsgMemoryGrp(keyval)
{
  if (debug_) cout << "MTMPIMemoryGrp keyval CTOR entered" << endl;

  th_ = ThreadGrp::get_default_threadgrp();

  KeyValValueint nthreaddef(th_->nthread() + 1);
  int nthread = keyval->intvalue("num_threads",nthreaddef);
  cout << node0 << indent << "MTMPIMemoryGrp: num_threads = " << nthread << endl;

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

  th_ = th_->clone(nthread);
  nthread = th_->nthread();
  if (nthread < 2) {
      cout << "MTMPIMemoryGrp didn't get enough threads" << endl;
      abort();
    }

  if (debug_) {
      char name[256];
      sprintf(name, "mpqc.hand.%d", me());
      hout.open(name);
      sprintf(name, "mpqc.main.%d", me());
      mout.open(name);
    }

  req_type_ = 9125;
  to_type_ = 9126;
  fr_type_ = 9127;

  thread_ = new Thread*[nthread-1];
  th_->add_thread(0,0);
  for (i=1; i<nthread; i++) {
      thread_[i-1] = new MTMPIThread(this,req_type_,to_type_,fr_type_);
      th_->add_thread(i,thread_[i-1]);
    }
  mem_lock_ = th_->new_lock();
  print_lock_ = th_->new_lock();
}

long
MTMPIMemoryGrp::lockcomm()
{
  mem_lock_->lock();
  return 0;
}

void
MTMPIMemoryGrp::unlockcomm(long oldvalue)
{
  mem_lock_->unlock();
}

void
MTMPIMemoryGrp::retrieve_data(void *data, int node, int offset, int size)
{
  MemoryDataRequest req(MemoryDataRequest::Retrieve,me(),offset,size);

  // send the request
  if (debug_) {
      print_lock_->lock();
      req.print("SEND",mout);
      print_lock_->unlock();
    }
  MPI_Send(req.data(),req.nbytes(),MPI_BYTE,node,req_type_,MPI_COMM_WORLD);

  // receive the data
  MPI_Status status;
  MPI_Recv(data,size,MPI_BYTE,node,fr_type_,MPI_COMM_WORLD,&status);
}

void
MTMPIMemoryGrp::replace_data(void *data, int node, int offset, int size)
{
  MemoryDataRequest req(MemoryDataRequest::Replace,me(),offset,size);

  if (debug_) {
      print_lock_->lock();
      req.print("SEND",mout);
      print_lock_->unlock();
    }
  MPI_Send(req.data(),req.nbytes(),MPI_BYTE,node,to_type_,MPI_COMM_WORLD);

  // wait for the go ahead message
  MPI_Status status;
  char junk;
  MPI_Recv(&junk,1,MPI_BYTE,node,fr_type_,MPI_COMM_WORLD,&status);

  // send the data
  MPI_Send(data,size,MPI_BYTE,node,to_type_,MPI_COMM_WORLD);

  // wait for the ack message
  MPI_Recv(&junk,1,MPI_BYTE,node,fr_type_,MPI_COMM_WORLD,&status);
}

void
MTMPIMemoryGrp::sum_data(double *data, int node, int offset, int size)
{
  MemoryDataRequest req(MemoryDataRequest::DoubleSum,me(),offset,size);

  // send the request
  if (debug_) {
      print_lock_->lock();
      req.print("SEND",mout);
      print_lock_->unlock();
    }
  MPI_Send(req.data(),req.nbytes(),MPI_BYTE,node,req_type_,MPI_COMM_WORLD);

  // wait for the go ahead message
  MPI_Status status;
  char junk;
  MPI_Recv(&junk,1,MPI_BYTE,node,fr_type_,MPI_COMM_WORLD,&status);

  int dsize = size/sizeof(double);
  int dremain = dsize;
  int dcurrent = 0;
  while(dremain>0) {
      int dchunksize = dbufsize;
      if (dremain < dchunksize) dchunksize = dremain;
      // send the data
      MPI_Send(&data[dcurrent],dchunksize,MPI_DOUBLE,
               node,to_type_,MPI_COMM_WORLD);
      dcurrent += dchunksize;
      dremain -= dchunksize;
    }

  // wait for the ack message
  MPI_Recv(&junk,1,MPI_BYTE,node,fr_type_,MPI_COMM_WORLD,&status);
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
      MPI_Send(req.data(),req.nbytes(),MPI_BYTE,me(),req_type_,MPI_COMM_WORLD);
    }

  // wait on the thread to shutdown
  th_->wait_threads();
}

#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
