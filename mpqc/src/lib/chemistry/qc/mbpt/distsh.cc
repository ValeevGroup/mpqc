//
// distsh.cc
// based on: csgrade12.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Ida Nielsen <ida@kemi.aau.dk>
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
#pragma implementation
#endif

#include <util/misc/formio.h>
#include <chemistry/qc/mbpt/distsh.h>

using namespace std;

/////////////////////////////////////////////////////////////////
// Function iquicksort performs a quick sort (larger -> smaller) 
// of the integer data in item by the integer indices in index;
// data in item remain unchanged
/////////////////////////////////////////////////////////////////
static void
iqs(int *item,int *index,int left,int right)
{
  register int i,j;
  int x,y;
 
  i=left; j=right;
  x=item[index[(left+right)/2]];
 
  do {
    while(item[index[i]]>x && i<right) i++;
    while(x>item[index[j]] && j>left) j--;
 
    if (i<=j) {
      if (item[index[i]] != item[index[j]]) {
        y=index[i];
        index[i]=index[j];
        index[j]=y;
        }
      i++; j--;
      }
    } while(i<=j);
       
  if (left<j) iqs(item,index,left,j);
  if (i<right) iqs(item,index,i,right);
}

static void
iquicksort(int *item,int *index,int n)
{
  int i;
  if (n<=0) return;
  for (i=0; i<n; i++) {
    index[i] = i;
    }
  iqs(item,index,0,n-1);
  }

/////////////////////////////////////////////////////////////////
// DistShellPair class

DistShellPair::DistShellPair(const Ref<MessageGrp> & msg,
                             int nthread, int mythread,
                             const Ref<ThreadLock> & lock,
                             const Ref<GaussianBasisSet> & basis):
  msg_(msg),
  nthread_(nthread),
  mythread_(mythread),
  lock_(lock),
  basis_(basis)
{
  ncpu_ = nthread_*msg->n();
  ncpu_less_0_ = nthread_*(msg->n()-1);

  print_percent_ = 10;
  debug_ = 0;
  dynamic_ = 0;
  req_type_ = 18101;
  ans_type_ = 18102;

  init();
}

DistShellPair::~DistShellPair()
{
}

void
DistShellPair::set_dynamic(int d)
{
  dynamic_ = d;
  if (msg_->n() <= 1) {
    dynamic_ = 0;
    }
}

void
DistShellPair::init()
{
  S_ = 0;
  R_ = msg_->me()*nthread_ + mythread_;
  while (R_ > S_) {
    S_++;
    R_ = R_ - S_;
    }

  int nshell = basis_->nshell();
  int ntri = (nshell*(nshell+1))/2;
  ntask_ = ntri/ncpu_;
  print_index_ = 0;
  print_interval_ = print_percent_*ntask_/100;
  if (print_interval_==0) print_interval_ = 1;
}

void
DistShellPair::serve_tasks()
{
  // intialize work arrays
  int S,R,index;
  int nshell = basis_->nshell();
  int ntri = (nshell*(nshell+1))/2;
  int *cost = new int[ntri];
  int *Svec = new int[ntri];
  int *Rvec = new int[ntri];
  int *Ivec = new int[ntri];
  index = 0;
  for (S=0; S<nshell; S++) {
    for (R=0; R<=S; R++) {
      cost[index] = basis_->shell(S).nfunction()*basis_->shell(R).nfunction();
      Svec[index] = S;
      Rvec[index] = R;
      Ivec[index] = index;
      index++;
      }
    }

  // sort work
  iquicksort(cost, Ivec, ntri);
  if (debug_ > 1) {
    ExEnv::out() << "costs of shell pairs" << endl;
    for (index=0; index<ntri; index++) {
      ExEnv::out() << scprintf(" (%d %d):%d",Svec[Ivec[index]],Rvec[Ivec[index]],
                       cost[Ivec[index]])
           << endl;
      }
    }

  // process requests
  int nreq = ntri + ncpu_less_0_;
  int iwork = 0;
  int print_index = 0;
  int print_interval = print_percent_*nreq/100;
  if (print_interval==0) print_interval = 1;
  int nreq_left = nreq;
  while (nreq_left) {
    int node;
    msg_->recvt(req_type_,&node,1);
    int SR[2];
    if (iwork < ntri) {
      SR[0] = Svec[Ivec[iwork]];
      SR[1] = Rvec[Ivec[iwork]];
      iwork++;
      if (print_index++%print_interval == 0) {
        ExEnv::out() << indent
             << scprintf("sending shell pair (%3d %3d) to %3d, %4.1f%% complete",
                         SR[0],SR[1],node,(double)(print_index*100)/nreq)
             << endl;
        }
      }
    else {
      SR[0] = -1;
      SR[1] = -1;
      if (print_index++%print_interval == 0) {
        ExEnv::out() << indent
             << scprintf("sending no more tasks message to %3d, %4.1f%% complete",
                         node,(double)(print_index*100)/nreq)
             << endl;
        }
      }
    msg_->sendt(node,ans_type_,SR,2);
    nreq_left--;
    }

  if (debug_) {
      ExEnv::out() << "all requests processed" << endl;
    }

  delete[] cost;
  delete[] Svec;
  delete[] Rvec;
  delete[] Ivec;
}

int
DistShellPair::get_task(int &S, int &R)
{
  if (dynamic_) { // dynamic load balancing
    int me = msg_->me();
    if (me == 0) {
      if (mythread_ == 0) serve_tasks();
      return 0;
      }
    else {
      int SR[2];

      lock_->lock();
      msg_->sendt(0,req_type_,&me,1);
      msg_->recvt(ans_type_,SR,2);
      lock_->unlock();

      S = SR[0];
      R = SR[1];
      if (S == -1) return 0;
      }
    }
  else { // static load balancing
    if (S_ >= basis_->nshell()) return 0;
    S = S_;
    R = R_;
    // advance to the next S_, R_
    R_ += ncpu_;
    while (R_ > S_) {
      S_++;
      R_ = R_ - S_;
      }
    if (print_index_++%print_interval_ == 0) {
      if (mythread_ == 0 && msg_->me() == 0) {
        ExEnv::out() << indent
             << scprintf("working on shell pair (%3d %3d), %4.1f%% complete",
                         S,R,(double)(print_index_*100)/ntask_)
             << endl;
        }
      }
    }
  return 1;
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
