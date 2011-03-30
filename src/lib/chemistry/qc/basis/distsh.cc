//
// distsh.cc
// modeled after distsh.cc
//
// Copyright (C) 2009 Edward F. Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
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

#include <math.h>
#include <algorithm>
#include <util/misc/formio.h>
#include <chemistry/qc/basis/distsh.h>

using namespace std;
using namespace sc;

// Defining REVERSE_ORDERING does the small tasks first.  This would give
// poorer load balancing and it included only for debugging (in particular,
// for stressing MPI libraries up front instead of at the end of a run).
#undef REVERSE_ORDERING

struct IndexPlusCost {
  IndexPlusCost() {}
  IndexPlusCost(int I, int C) : i(I), cost(C) {}
  int i;
  int cost;
};
struct DecreasingCost {
  bool operator()(const IndexPlusCost& A,
                  const IndexPlusCost& B) {
    return A.cost > B.cost;
  }
};
struct IncreasingCost {
  bool operator()(const IndexPlusCost& A,
                  const IndexPlusCost& B) {
    return A.cost > B.cost;
  }
};

/////////////////////////////////////////////////////////////////
// DistShell class

DistShell::DistShell(const Ref<MessageGrp> & msg,
                     int nthread, int mythread,
                     const Ref<ThreadLock> & lock,
                     const Ref<GaussianBasisSet> & bs,
                     bool dynamic,
                     int max_nfunction,
                     int max_nshell,
                     SharedData *shared):
  msg_(msg),
  nthread_(nthread),
  mythread_(mythread),
  lock_(lock),
  bs_(bs),
  max_nfunction_(max_nfunction),
  max_nshell_(max_nshell),
  task_dynamic_(dynamic),
  thread_dynamic_(false),
  shared_(shared)
{
  // validate input
  if (max_nfunction < bs->max_nfunction_in_shell()) {
    max_nfunction = bs->max_nfunction_in_shell()*max_nshell;
  }
  if (max_nshell < 1) max_nshell = 1;

  ncpu_ = nthread_*msg->n();

  // Cannot do dynamic load balancing if there's only 1 task
  if (msg->n() == 1) {
    task_dynamic_ = false;
    if (dynamic && shared_ != 0) thread_dynamic_ = true;
  }
  debug_ = 0;
  print_percent_ = 10;
  req_type_ = 18101;
  ans_type_ = 18102;

  init_work();
  set_print_percent(print_percent_);

  init();
}

DistShell::~DistShell()
{
}

void
DistShell::init()
{
  current_task_ = 0;
}

void
DistShell::init_work()
{
  const long int nsh = bs_->nshell();

  // sort shells into decreasing shell sizes
  {
    std::vector<IndexPlusCost> shell_indices(nsh);
    for(int s=0; s<nsh; ++s) {
      shell_indices[s] = IndexPlusCost(s, bs_->shell(s).nfunction());
    }
    std::stable_sort(shell_indices.begin(), shell_indices.end(), DecreasingCost());
    shell_indices.resize(nsh);
    for(int s=0; s<nsh; ++s)
      shell_map_[s] = shell_indices[s].i;
  }

  // divide shells into work units
  {
    int shell = 0;
    int current_payload_nfunctions = 0;
    int current_payload_nshells = 0;
    while(shell < nsh) {
      Task current_task;  current_task.first = shell;
      while (current_payload_nfunctions <= max_nfunction_ &&
             current_payload_nshells <= max_nshell_) {
        current_payload_nfunctions += bs_->shell(shell_map_[shell]).nfunction();
        ++current_payload_nshells;
        ++shell;
      }
      current_task.second = current_payload_nshells;
      tasks_.push_back(current_task);
      current_payload_nfunctions = 0;
      current_payload_nshells = 0;
    }
  }

  ntask_ = tasks_.size();
}

void
DistShell::serve_tasks()
{
  // process requests
  const long int ncpu_less_0 = nthread_*(msg_->n() - 1);
  long int nreq = ntask_ + ncpu_less_0;
  int nreq_left = nreq;
  while (nreq_left) {
    int node;
    msg_->recvt(MessageGrp::AnySender,req_type_,&node,1);
    if (current_task_ < ntask_) {
      msg_->sendt(node,ans_type_,&current_task_,1);
      if (print_percent_ <= 100.0
          && current_task_%print_interval_ == 0) {
        ExEnv::outn() << indent
                      << scprintf("sent shell (%3d) to %3d, %6.3f%% complete",
                                  shell_map_[tasks_[current_task_].first],
                                  node,
                                  (double)current_task_*100.0/nreq)
                      << " (" << current_task_ << " of " << ntask_ << ")"
                      << endl;
      }
      current_task_++;
    }
    else {
      const int end_of_tasks = -1;
      msg_->sendt(node,ans_type_,&end_of_tasks,1);
      if (print_percent_ <= 100.0
          && current_task_%print_interval_ == 0) {
        ExEnv::outn() << indent
		     << scprintf("sent \"no more tasks\" message to %3d, %6.3f%% complete",
				 node,(double)current_task_*100.0/nreq)
		     << endl;
      }
      current_task_++;
    }
    nreq_left--;
  }

  if (debug_) {
    ExEnv::outn() << "all requests processed" << endl;
  }
}

int
DistShell::get_task(int &N, int &I)
{
  const int me = msg_->me();
  if (task_dynamic_) { // dynamic load balancing
    if (me == 0) {
      if (mythread_ == 0) serve_tasks();
      return 0;
    }
    else {
      int n;

      lock_->lock();
      msg_->sendt(0,req_type_,&me,1);
      msg_->recvt(MessageGrp::AnySender,ans_type_,&n,1);
      lock_->unlock();

      if (n == -1) return 0;
      N = tasks_[n].first;
      I = tasks_[n].second;
    }
  }
  else if (thread_dynamic_) {
    lock_->lock();
    const int my_shell = shared_->shell_++;
    lock_->unlock();
    if (my_shell < ntask_) {
      N = tasks_[my_shell].first;
      I = tasks_[my_shell].second;
    }
    else {
      return 0;
    }
    if (my_shell%print_interval_ == 0) {
      if (print_percent_ <= 100.0 && msg_->me() == 0) {
        ExEnv::outn() << indent
             << scprintf("  working on shell (%3d), %6.3f%% complete",
                 shell_map_[N],((double)my_shell*100.0)/ntask_)
                     << " (" << my_shell << " of " << ntask_ << ")"
             << endl;
      }
    }
  }
  else { // static load balancing
    while (current_task_%ncpu_ != me) ++current_task_;
    if (current_task_ < ntask_) {
      N = tasks_[current_task_].first;
      I = tasks_[current_task_].second;
    }
    else
      return 0;

    if (current_task_%print_interval_ == 0) {
      if (print_percent_ <= 100.0
          && mythread_ == 0 && msg_->me() == 0) {
        ExEnv::outn() << indent
		     << scprintf("  working on shell (%3d), %6.3f%% complete",
		                 shell_map_[N],((double)current_task_*100.0)/ntask_)
                     << " (" << current_task_ << " of " << ntask_ << ")"
		     << endl;
      }
    }
    current_task_++;
  }
  return 1;
}

void
DistShell::set_print_percent(double pi)
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
