//
// distsh.cc
// based on: mbpt/distsh.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Ida Nielsen <ida@kemi.aau.dk>
// Updated: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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

#include <math.h>

#include <util/misc/formio.h>
#include <chemistry/qc/basis/distshpair.h>

using namespace std;
using namespace sc;

// Defining REVERSE_ORDERING does the small tasks first.  This would give
// poorer load balancing and it included only for debugging (in particular,
// for stressing MPI libraries up front instead of at the end of a run).
#undef REVERSE_ORDERING

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
iquicksort(int *item,int *index,long int n)
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
			     const Ref<GaussianBasisSet> & bs1, const Ref<GaussianBasisSet> & bs2,
			     bool dynamic, SharedData *shared):
  msg_(msg),
  nthread_(nthread),
  mythread_(mythread),
  lock_(lock),
  bs1_(bs1), bs2_(bs2),
  task_dynamic_(dynamic),
  thread_dynamic_(false),
  cost_(0),
  Svec_(0),
  Rvec_(0),
  Ivec_(0),
  shared_(shared)
{
  ncpu_ = nthread_*msg->n();
  ncpu_less_0_ = nthread_*(msg->n()-1);

  // Cannot do dynamic load balancing if there's only 1 task
  if (msg->n() == 1) {
    task_dynamic_ = false;
    if (dynamic && shared_ != 0) thread_dynamic_ = true;
    }
  debug_ = 0;
  print_percent_ = 10;
  bs1_eq_bs2_ = (bs1_ == bs2_);
  req_type_ = 18101;
  ans_type_ = 18102;

  // Only for static case with different basis sets
  if (!bs1_eq_bs2_) {
    int nsh2 = bs2_->nshell();
    incS_ = ncpu_/nsh2;
    incR_ = ncpu_%nsh2;
  }

  long int nsh1 = bs1_->nshell();
  long int nsh2 = bs2_->nshell();
  long int nshpairs = (bs1_eq_bs2_ ? (nsh1*(nsh1+1))/2 : nsh1*nsh2);
  if (task_dynamic_ || thread_dynamic_) {
    ntask_ = nshpairs;
    init_dynamic_work();
    }
  else {
    ntask_ = nshpairs/ncpu_;
    }
  set_print_percent(print_percent_);

  init();
}

DistShellPair::~DistShellPair()
{
  delete[] cost_;
  delete[] Svec_;
  delete[] Rvec_;
  delete[] Ivec_;
}

void
DistShellPair::init()
{
    // Compute starting S_ and R_ for this thread
    if (bs1_eq_bs2_) {
      // This is a slightly nonobvious computation of S_ and R_ from SR = msg_->me()*nthread_ + mythread_
      // when bs1_eq_bs2 == true
      S_ = 0;
      R_ = msg_->me()*nthread_ + mythread_;
      while (R_ > S_) {
        S_++;
        R_ = R_ - S_;
        }
      }
    else {
      // Things are simple when basis sets are different
      long int nsh2 = bs2_->nshell();
      long int absthreadindex = msg_->me()*nthread_ + mythread_;
      S_ = absthreadindex/nsh2;
      R_ = absthreadindex%nsh2;
      }
  current_shellpair_ = 0;
}

void
DistShellPair::init_dynamic_work()
{
  // initialize work arrays
  int S,R,index;
  int nsh1 = bs1_->nshell();
  int nsh2 = bs2_->nshell();
  delete[] cost_;
  delete[] Svec_;
  delete[] Rvec_;
  delete[] Ivec_;
  cost_ = new int[ntask_];
  Svec_ = new int[ntask_];
  Rvec_ = new int[ntask_];
  Ivec_ = new int[ntask_];
  index = 0;
  for (S=0; S<nsh1; S++) {
    int Rmax = (bs1_eq_bs2_ ? S : nsh2-1);
    for (R=0; R<=Rmax; R++) {
      cost_[index] = bs1_->shell(S).nfunction() * bs2_->shell(R).nfunction() *
	            bs1_->shell(S).nprimitive() * bs2_->shell(R).nprimitive();
#ifdef REVERSE_ORDERING
      cost_[index] = - cost[index];
#endif
      Svec_[index] = S;
      Rvec_[index] = R;
      Ivec_[index] = index;
      index++;
    }
  }

  // sort work
  iquicksort(cost_, Ivec_, ntask_);
  if (debug_ > 1) {
    ExEnv::outn() << "costs of shell pairs" << endl;
    for (index=0; index<ntask_; index++) {
      ExEnv::outn() << scprintf(" (%d %d):%d",Svec_[Ivec_[index]],Rvec_[Ivec_[index]],
			       cost_[Ivec_[index]])
		   << endl;
    }
  }
}

void
DistShellPair::serve_tasks()
{
  // process requests
  long int nreq = ntask_ + ncpu_less_0_;
  int nreq_left = nreq;
  while (nreq_left) {
    int node;
    msg_->recvt(MessageGrp::AnySender,req_type_,&node,1);
    int SR[2];
    if (current_shellpair_ < ntask_) {
      SR[0] = Svec_[Ivec_[current_shellpair_]];
      SR[1] = Rvec_[Ivec_[current_shellpair_]];
      msg_->sendt(node,ans_type_,SR,2);
      if (print_percent_ <= 100.0
          && current_shellpair_%print_interval_ == 0) {
        ExEnv::outn() << indent
                      << scprintf("sent shell pair (%3d %3d) to %3d, %6.3f%% complete",
                                  SR[0],SR[1],node,(double)current_shellpair_*100.0/nreq)
                      << " (" << current_shellpair_ << " of " << ntask_ << ")"
                      << endl;
      }
      current_shellpair_++;
    }
    else {
      SR[0] = -1;
      SR[1] = -1;
      msg_->sendt(node,ans_type_,SR,2);
      if (print_percent_ <= 100.0
          && current_shellpair_%print_interval_ == 0) {
        ExEnv::outn() << indent
		     << scprintf("sent no more tasks message to %3d, %6.3f%% complete",
				 node,(double)current_shellpair_*100.0/nreq)
		     << endl;
      }
      current_shellpair_++;
    }
    nreq_left--;
  }

  if (debug_) {
    ExEnv::outn() << "all requests processed" << endl;
  }
}

int
DistShellPair::get_task(int &S, int &R)
{
  if (task_dynamic_) { // dynamic load balancing
    int me = msg_->me();
    if (me == 0) {
      if (mythread_ == 0) serve_tasks();
      return 0;
    }
    else {
      int SR[2];
      
      lock_->lock();
      msg_->sendt(0,req_type_,&me,1);
      msg_->recvt(MessageGrp::AnySender,ans_type_,SR,2);
      lock_->unlock();

      S = SR[0];
      R = SR[1];
      if (S == -1) return 0;
    }
  }
  else if (thread_dynamic_) {
    long my_shellpair;
    lock_->lock();
    my_shellpair = shared_->shellpair_++;
    lock_->unlock();
    if (my_shellpair < ntask_) {
      S = Svec_[Ivec_[my_shellpair]];
      R = Rvec_[Ivec_[my_shellpair]];
      }
    else {
      S = R = -1;
      return 0;
      }
    if (my_shellpair%print_interval_ == 0) {
      if (print_percent_ <= 100.0 && msg_->me() == 0) {
        ExEnv::outn() << indent 
             << scprintf("  working on shell pair (%3d %3d), %6.3f%% complete",
                 S,R,((double)my_shellpair*100.0)/ntask_)
                     << " (" << my_shellpair << " of " << ntask_ << ")"
             << endl;
      }
    }
  }
  else { // static load balancing
    int nsh1 = bs1_->nshell();
    int nsh2 = bs2_->nshell();
    if (S_ >= nsh1 || R_ >= nsh2) return 0;
    S = S_;
    R = R_;
    // advance to the next S_, R_
    if (bs1_eq_bs2_) {
      R_ += ncpu_;
      while (R_ > S_) {
	S_++;
	R_ = R_ - S_;
      }
    }
    else {
      S_ += incS_;
      R_ += incR_;
      if (R_ >= nsh2) {	S_++;
	R_ -= nsh2;
      }
    }
    if (current_shellpair_%print_interval_ == 0) {
      if (print_percent_ <= 100.0
          && mythread_ == 0 && msg_->me() == 0) {
        ExEnv::outn() << indent 
		     << scprintf("  working on shell pair (%3d %3d), %6.3f%% complete",
				 S,R,((double)current_shellpair_*100.0)/ntask_)
                     << " (" << current_shellpair_ << " of " << ntask_ << ")"
		     << endl;
      }
    }
    current_shellpair_++;
  }
  return 1;
}

void
DistShellPair::set_print_percent(double pi)
{
  print_percent_ = pi;
  double print_interval = (double) print_percent_ * ntask_ / 100;
  print_interval_ = (long int)floor(print_interval + 0.5);
  if (print_interval_ < 1)
    print_interval_ = 1;
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
