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

#include <cassert>

#include <util/misc/formio.h>
#include <util/group/messmpi.h>
#include <util/group/memmtmpi.h>

#define OMPI_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#include <mpi.h>

using namespace std;

// Define this to use immediate mode.  This was added added to work
// around bugs in non-immediate mode optimizations in an MPI impl.
#undef USE_IMMEDIATE_MODE

namespace sc {

static const int dbufsize = 32768;

///////////////////////////////////////////////////////////////////////
// The MTMPIThread class

class MTMPIThread: public Thread {
  private:
    MTMPIMemoryGrp *mem_;
    int req_tag_;
    int tag_;
    unsigned int nreq_recd_;
    double chunk[dbufsize];
    Ref<RegionTimer> timer_;
  public:
    MTMPIThread(MTMPIMemoryGrp *, int reqtype, int tag, bool use_timer);
    void run();
    int run_one();
    unsigned int nreq_recd() { return nreq_recd_; }
    void set_nreq_recd(unsigned int val) { nreq_recd_ = val; }
    Ref<RegionTimer> timer() { return timer_; }
};

MTMPIThread::MTMPIThread(MTMPIMemoryGrp *mem,
                         int reqtype, int tag,
                         bool use_timer)
{
  mem_ = mem;
  req_tag_ = reqtype;
  tag_ = tag;
  nreq_recd_ = 0;
  if (use_timer) timer_ = new RegionTimer("MTMPIThread",1,1);
}

void
MTMPIThread::run()
{
  while(run_one());
}

int
MTMPIThread::run_one()
{
  int i;
  long dsize;
  long dremain;
  long doffset;
  long l;
  MemoryDataRequest req;
  MPI_Status status;
  Timer tim(timer_,"recv_req");
#ifndef USE_IMMEDIATE_MODE
  MPI_Recv(req.data(),req.nbytes(),MPI_BYTE,MPI_ANY_SOURCE,
           req_tag_,mem_->comm_comm_,&status);
#else
  MPI_Request mpireq;
  MPI_Irecv(req.data(),req.nbytes(),MPI_BYTE,MPI_ANY_SOURCE,
           req_tag_,mem_->comm_comm_,&mpireq);
  MPI_Wait(&mpireq,&status);
#endif // USE_IMMEDIATE_MODE
  tim.change("handle_req");
  int rtag = req.serial_number();
  if (mem_->debug_) {
      mem_->print_lock_->lock();
      req.print("RECV",mem_->hout);
      mem_->print_lock_->unlock();
    }
  if (req.touches_data()) {
      MPQC_ASSERT(req.size() >= 0);
      if (req.offset()+req.size() > mem_->localsize()) {
          mem_->print_lock_->lock();
          req.print("BAD RECV");
          ExEnv::outn() << "mem_->localsize() = " << mem_->localsize() << endl;
          mem_->print_lock_->lock();
        }
      MPQC_ASSERT(req.offset()+req.size() <= mem_->localsize());
    }
  Timer tim2(timer_);
  if (timer_.nonnull()) tim2.enter(req.request_string());
  Timer tim3(timer_);
  switch (req.request()) {
  case MemoryDataRequest::Deactivate:
      return 0;
  case MemoryDataRequest::Retrieve:
      nreq_recd_++;
      if (req.lock())
          mem_->obtain_local_lock(req.offset(), req.offset()+req.size());
      MPI_Send(&mem_->data_[req.offset()],req.size(),MPI_BYTE,
               req.node(),rtag,mem_->comp_comm_);
      break;
  case MemoryDataRequest::Replace:
      nreq_recd_++;
      // May be able to get rid of this MPI_Send - MLL
      MPI_Send(&tag_,1,MPI_INT,req.node(),rtag,mem_->comp_comm_);
      MPI_Recv(&mem_->data_[req.offset()],req.size(),MPI_BYTE,
               req.node(),tag_,mem_->comm_comm_,&status);
      if (req.lock())
          mem_->release_local_lock(req.offset(), req.offset()+req.size());
      break;
  case MemoryDataRequest::DoubleSum:
      nreq_recd_++;
//      MPI_Send(&tag_,1,MPI_INT,req.node(),rtag,mem_->comm_);
      dsize = req.size()/sizeof(double);
      dremain = dsize;
      doffset = req.offset()/sizeof(double);
      tim3.enter("obtain_local_lock");
      mem_->obtain_local_lock(req.offset(), req.offset()+req.size());
      tim3.exit();
      while(dremain>0) {
          int dchunksize = dbufsize;
          if (dremain < dchunksize) dchunksize = dremain;
          tim3.enter("recv");
#ifndef USE_IMMEDIATE_MODE
          MPI_Recv(chunk,dchunksize,MPI_DOUBLE,
                   req.node(),rtag ,mem_->comm_comm_,&status);
#else
          MPI_Request mpireq;
          MPI_Irecv(chunk,dchunksize,MPI_DOUBLE,
                   req.node(),rtag ,mem_->comm_comm_,&mpireq);
          MPI_Wait(&mpireq,&status);
#endif // USE_IMMEDIATE_MODE
          tim3.change("sum");
          double *source_data = &((double*)mem_->data_)[doffset];
          for (i=0; i<dchunksize; i++) {
              source_data[i] += chunk[i];
            }
          tim3.exit();
          dremain -= dchunksize;
          doffset += dchunksize;
        }
      tim3.enter("release_local_lock");
      mem_->release_local_lock(req.offset(), req.offset()+req.size());
      tim3.exit();
      break;
  default:
      mem_->print_lock_->lock();
      ExEnv::outn() << "MTMPIThread: bad memory data request" << endl;
      mem_->print_lock_->unlock();
      abort();
    }
  return 1;
}

///////////////////////////////////////////////////////////////////////
// The MTMPIMemoryGrp class

static ClassDesc MTMPIMemoryGrp_cd(
  typeid(MTMPIMemoryGrp),"MTMPIMemoryGrp",1,"public ActiveMsgMemoryGrp",
  0, create<MTMPIMemoryGrp>, 0);

MTMPIMemoryGrp::MTMPIMemoryGrp(const Ref<MessageGrp>& msg,
                               const Ref<ThreadGrp>& th,
                               MPI_Comm comm):
  ActiveMsgMemoryGrp(msg),
  nbuffer_(0)
{
  if (debug_) ExEnv::outn() << "MTMPIMemoryGrp CTOR entered" << endl;

  th_ = th;

  init_mtmpimg(comm,th_->nthread());
}

MTMPIMemoryGrp::MTMPIMemoryGrp(const Ref<KeyVal>& keyval):
  ActiveMsgMemoryGrp(keyval)
{
  if (debug_) ExEnv::outn() << "MTMPIMemoryGrp keyval CTOR entered" << endl;

  th_ = ThreadGrp::get_default_threadgrp();

  KeyValValueint nthreaddef(th_->nthread());
  int nthread = keyval->intvalue("num_threads",nthreaddef);
  ExEnv::out0() << indent << "MTMPIMemoryGrp: num_threads = " << nthread << endl;

  KeyValValueint nbufferdef(0);
  nbuffer_ = keyval->intvalue("num_buffer",nbufferdef);
  ExEnv::out0() << indent << "MTMPIMemoryGrp: num_buffer = "
                << nbuffer_ << endl;

  KeyValValueboolean usetimerdef(false);
  bool use_timer = keyval->booleanvalue("use_timer",usetimerdef);
  ExEnv::out0() << indent << "MTMPIMemoryGrp: use_timer = "
                << use_timer << endl;
  if (use_timer) {
      timer_ = new RegionTimer("MTMPIMemoryGrp",1,1);
    }

  init_mtmpimg(MPI_COMM_WORLD,nthread);
}

MTMPIMemoryGrp::~MTMPIMemoryGrp()
{
  deactivate();
  for (int i=0; i<th_->nthread()-1; i++) {
      delete thread_[i];
    }
  delete[] thread_;
  delete[] nreq_sent_;
  delete[] nreq_sent_buf_;
}

void
MTMPIMemoryGrp::init_mtmpimg(MPI_Comm comm, int nthread)
{
  int i;
  active_ = 0;

  if (nthread < 2) nthread = 2;
  th_ = th_->clone(nthread);
  nthread = th_->nthread();

  if (debug_) {
      char name[256];
      sprintf(name, "mpqc.hand.%d", me());
      hout.open(name);
      sprintf(name, "mpqc.main.%d", me());
      mout.open(name);
    }

  MPI_Comm_dup(comm, &comp_comm_);
  MPI_Comm_dup(comm, &comm_comm_);

  MPI_Errhandler_set(comp_comm_, MPI_ERRORS_ARE_FATAL);
  MPI_Errhandler_set(comm_comm_, MPI_ERRORS_ARE_FATAL);

  serial_ = 0;
  req_tag_ = 15001;

  serial_lock_ = th_->new_lock();

  thread_ = new MTMPIThread*[nthread-1];
  th_->add_thread(0,0);
  for (i=1; i<nthread; i++) {
      thread_[i-1] = new MTMPIThread(this,req_tag_,req_tag_ + i,
                                     timer_.nonnull());
      th_->add_thread(i,thread_[i-1]);
    }
  print_lock_ = th_->new_lock();

  nreq_sent_ = new unsigned int[n()];
  memset(nreq_sent_, 0, sizeof(unsigned int)*n());
  nreq_sent_buf_ = new unsigned int[n()];
}

int
MTMPIMemoryGrp::serial(int node)
{
  serial_lock_->lock();
  nreq_sent_[node]++;
  int r = serial_;
  serial_++;
  if (serial_ == req_tag_) serial_ = 0;
  serial_lock_->unlock();
  return r;
}

void
MTMPIMemoryGrp::retrieve_data(void *data, int node, long offset, long size,
                              int lock)
{
  MemoryDataRequest req(MemoryDataRequest::Retrieve,me(),offset,size,lock,
                        serial(node));
  int tag = req.serial_number();

  // send the request
  if (debug_) {
      print_lock_->lock();
      req.print("SEND",mout);
      print_lock_->unlock();
    }
  MPI_Send(req.data(),req.nbytes(),MPI_BYTE,node,req_tag_,comm_comm_);

  // receive the data
  MPI_Status status;
  MPI_Recv(data,size,MPI_BYTE,node,tag,comp_comm_,&status);
}

void
MTMPIMemoryGrp::replace_data(void *data, int node, long offset, long size,
                             int unlock)
{
  MemoryDataRequest req(MemoryDataRequest::Replace,me(),offset,size,unlock,
                        serial(node));
  int tag = req.serial_number();

  if (debug_) {
      print_lock_->lock();
      req.print("SEND",mout);
      print_lock_->unlock();
    }
  MPI_Send(req.data(),req.nbytes(),MPI_BYTE,node,req_tag_,comm_comm_);

  // wait for the go ahead message
  MPI_Status status;
  int rtag;
  MPI_Recv(&rtag,1,MPI_INT,node,tag,comp_comm_,&status);

  // send the data
  MPI_Send(data,size,MPI_BYTE,node,rtag,comm_comm_);
}

void
MTMPIMemoryGrp::sum_data(double *data, int node, long offset, long size)
{
  MemoryDataRequest req(MemoryDataRequest::DoubleSum,me(),offset,size,
                        0, serial(node));
  int tag = req.serial_number();

  // send the request
  if (debug_) {
      print_lock_->lock();
      req.print("SEND",mout);
      print_lock_->unlock();
    }

  if (nbuffer_ > 0) {
      ThreadLockHolder lock(buffer_lock_);
      Timer tim1(timer_, "sum_data(request)");
      Timer tim2(timer_, "get_request");
      int bufnum = get_request();
      tim2.change("memcpy");
      memcpy(datareqs_[bufnum].data(),req.data(),req.nbytes());
      tim2.change("isend");
      MPI_Isend(datareqs_[bufnum].data(),req.nbytes(),MPI_BYTE,
                node,req_tag_,comm_comm_,&datareqs_mpireq_[bufnum]);
    }
  else {
#ifndef USE_IMMEDIATE_MODE
      MPI_Send(req.data(),req.nbytes(),MPI_BYTE,node,req_tag_,comm_comm_);
#else
      MPI_Status status;
      MPI_Request mpireq;
      MPI_Isend(req.data(),req.nbytes(),MPI_BYTE,
                node,req_tag_,comm_comm_,&mpireq);
      MPI_Wait(&mpireq,&status);
#endif // USE_IMMEDIATE_MODE
    }

  // wait for the go ahead message
//  MPI_Status status;
//  int rtag;
//  MPI_Recv(&rtag,1,MPI_INT,node,tag,comm_,&status);

  int dsize = size/sizeof(double);
  int dremain = dsize;
  int dcurrent = 0;
  while(dremain>0) {
      int dchunksize = dbufsize;
      if (dremain < dchunksize) dchunksize = dremain;
      // send the data
      if (nbuffer_ > 0) {
          ThreadLockHolder lock(buffer_lock_);
          Timer tim1(timer_, "sum_data(data)");
          Timer tim2(timer_, "get_buffer");
          int bufnum = get_buffer();
          tim2.change("memcpy");
          memcpy(databufs_[bufnum],&data[dcurrent],dchunksize*sizeof(double));
          tim2.change("isend");
          MPI_Isend(databufs_[bufnum],dchunksize,MPI_DOUBLE,
                    node,tag,comm_comm_,&databufs_mpireq_[bufnum]);
        }
      else {
#ifndef USE_IMMEDIATE_MODE
          MPI_Send(&data[dcurrent],dchunksize,MPI_DOUBLE,
                   node,tag,comm_comm_);
#else
          MPI_Request mpireq;
          MPI_Isend(&data[dcurrent],dchunksize,MPI_DOUBLE,
                    node,tag,comm_comm_,&mpireq);
          MPI_Wait(&mpireq,&status);
#endif // USE_IMMEDIATE_MODE
        }
      dcurrent += dchunksize;
      dremain -= dchunksize;
    }
}

void
MTMPIMemoryGrp::activate()
{
  // Only remote requests require the handler.  There are only remote
  // requests if there is more than one node.
  if (n() == 1) return;

  if (th_->nthread() < 2) {
      ExEnv::outn() << "MTMPIMemoryGrp didn't get enough threads" << endl;
      abort();
    }

  if (active_) return;
  active_ = 1;

  if (timer_.nonnull()) {
      timer_->reset();
      for (int i=0; i<th_->nthread()-1; i++) {
          thread_[i]->timer()->reset();
        }
    }

  if (nbuffer_ > 0) init_buffer();

  th_->start_threads();
}

void
MTMPIMemoryGrp::deactivate()
{
  if (!active_) return;
  active_ = 0;

  if (nbuffer_ > 0) {
      done_buffers();
    }

  // send a shutdown message
  MemoryDataRequest req(MemoryDataRequest::Deactivate);
  if (debug_) {
      print_lock_->lock();
      req.print("SEND",mout);
      print_lock_->unlock();
    }
  for (int i=1; i<th_->nthread(); i++) {
#ifndef USE_IMMEDIATE_MODE
      MPI_Send(req.data(),req.nbytes(),MPI_BYTE,me(),req_tag_,comm_comm_);
#else
      MPI_Request mpireq;
      MPI_Status status;
      MPI_Isend(req.data(),req.nbytes(),MPI_BYTE,me(),req_tag_,comm_comm_,&mpireq);
      MPI_Wait(&mpireq,&status);
#endif // USE_IMMEDIATE_MODE
    }

  // wait on the thread to shutdown
  th_->wait_threads();

  if (timer_.nonnull()) {
      for (int i=0; i<th_->nthread()-1; i++) {
          timer_->merge(thread_[i]->timer());
        }
      Ref<RegionTimer> global_timer = RegionTimer::default_regiontimer();
      if (global_timer.nonnull()) {
          global_timer->merge(timer_);
        }
    }
}

void
MTMPIMemoryGrp::sync()
{
  if (active_) {
      MPI_Allreduce(nreq_sent_, nreq_sent_buf_,
                    n(), MPI_UNSIGNED, MPI_SUM, comm_comm_);
      deactivate();
      unsigned int nreq_recd = 0;
      for (int i=0; i<th_->nthread()-1; i++) {
          nreq_recd += thread_[i]->nreq_recd();
          thread_[i]->set_nreq_recd(0);
        }
      int n_outstanding = nreq_sent_buf_[me()] - nreq_recd;
      for (int i=0; i<n_outstanding; i++) {
          thread_[0]->run_one();
        }
      memset(nreq_sent_, 0, sizeof(unsigned int)*n());
      // Make sure processing of all outstanding requests is finished
      // before starting the next phase.
      MPI_Barrier(comm_comm_);
      activate();
    }
  else {
      MPI_Barrier(comm_comm_);
    }
}

// nbuffer_ must be set and > 0 before this is called.
void
MTMPIMemoryGrp::init_buffer()
{
  current_datareq_index_ = 0;
  datareqs_.resize(nbuffer_);
  datareqs_mpireq_.resize(nbuffer_);
  std::fill(datareqs_mpireq_.begin(), datareqs_mpireq_.end(),
            MPI_REQUEST_NULL);

  current_data_index_ = 0;
  databufs_.resize(nbuffer_);
  databufs_mpireq_.resize(nbuffer_);
  std::fill(databufs_mpireq_.begin(), databufs_mpireq_.end(),
            MPI_REQUEST_NULL);

  for (int i=0; i<nbuffer_; i++) databufs_[i] = new double[dbufsize];

  buffer_lock_ = th_->new_lock();
}

// The lock must be held from before this is called until
// after the immediate mode MPI call using the buffer is made.
int
MTMPIMemoryGrp::next_buffer(int &counter,
                            std::vector<MPI_Request> &reqs)
{
  int rc = counter;
  counter++;
  if (counter>=nbuffer_) counter=0;

  if (reqs[rc] != MPI_REQUEST_NULL) {
      MPI_Status stat;
      MPI_Wait(&reqs[rc], &stat);
      reqs[rc] = MPI_REQUEST_NULL;
    }
  
  return rc;
}

int
MTMPIMemoryGrp::get_buffer()
{
  int bufnum = next_buffer(current_data_index_,
                           databufs_mpireq_);
  return bufnum;
}

int
MTMPIMemoryGrp::get_request()
{
  int bufnum = next_buffer(current_datareq_index_,
                           datareqs_mpireq_);
  return bufnum;
}

void
MTMPIMemoryGrp::done_buffers()
{
  std::vector<MPI_Status> stats(nbuffer_);
  MPI_Waitall(nbuffer_, &databufs_mpireq_.front(), &stats.front());
  MPI_Waitall(nbuffer_, &datareqs_mpireq_.front(), &stats.front());

  datareqs_.resize(0);
  datareqs_mpireq_.resize(0);

  for (int i=0; i<nbuffer_; i++) delete[] databufs_[i];

  databufs_.resize(0);
  databufs_mpireq_.resize(0);
}

Ref<MemoryGrp>
MTMPIMemoryGrp::clone()
{
  if (class_desc() != ClassDesc::name_to_class_desc("MTMPIMemoryGrp")) {
      // this will throw
      return MemoryGrp::clone();
    }

  Ref<MemoryGrp> ret;
  ret = new MTMPIMemoryGrp(msg_->clone(), th_->clone(), comp_comm_);

  return ret;
}

#endif

/////////////////////////////////////////////////////////////////////////////

}

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
