//
// thpthd.cc
//
// Copyright (C) 1997 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
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

#ifdef HAVE_CONFIG_H
#include <mpqc_config.h>
#endif

#if HAVE_PTHREAD_H
#include <pthread.h>
#endif

#include <string.h>

#include <util/keyval/keyval.h>
#include <util/group/thpthd.h>
#include <util/misc/formio.h>

using namespace std;
using namespace sc;

namespace sc {

/////////////////////////////////////////////////////////////////////////////
// PthreadThreadLock class

class PthreadThreadLock : public ThreadLock {
  private:
    pthread_mutex_t mutex_;
    pthread_mutexattr_t attr_;
    
  public:
    PthreadThreadLock() {
      pthread_mutexattr_init(&attr_);
//#if defined(PTHREAD_MUTEX_FAST_NP)
//      pthread_mutexattr_setkind_np(&attr_, PTHREAD_MUTEX_FAST_NP);
//#elif defined(MUTEX_FAST_NP)
//      pthread_mutexattr_setkind_np(&attr_, MUTEX_FAST_NP);
//#endif
      pthread_mutex_init(&mutex_, &attr_);
    }

    ~PthreadThreadLock() {
      pthread_mutexattr_destroy(&attr_);
      pthread_mutex_destroy(&mutex_);
    }

    void lock() { pthread_mutex_lock(&mutex_); }
    void unlock() { pthread_mutex_unlock(&mutex_); }
};

/////////////////////////////////////////////////////////////////////////////
// PthreadThreadGrp members

static ClassDesc PthreadThreadGrp_cd(
  typeid(PthreadThreadGrp),"PthreadThreadGrp",1,"public ThreadGrp",
  0, create<PthreadThreadGrp>, 0);

PthreadThreadGrp::PthreadThreadGrp()
  : ThreadGrp()
{
  pthreads_ = new pthread_t[nthread_];
  init_attr();
}


PthreadThreadGrp::PthreadThreadGrp(const PthreadThreadGrp &tg,int nthread):
  ThreadGrp(tg, nthread)
{
  pthreads_ = new pthread_t[nthread_];
  init_attr();
}

PthreadThreadGrp::PthreadThreadGrp(const Ref<KeyVal>& keyval)
  : ThreadGrp(keyval)
{
  pthreads_ = new pthread_t[nthread_];
  init_attr();
}

PthreadThreadGrp::~PthreadThreadGrp()
{
  if (pthreads_) {
    delete[] pthreads_;
    pthreads_ = 0;
    delete[] attr_;
 }
//  delete attr_;
}

void
PthreadThreadGrp::init_attr()
{
  attr_ = new pthread_attr_t[nthread_];

  for (int i=0; i<nthread_; i++) {
    pthread_attr_init(&attr_[i]);
#if defined(PTHREAD_CREATE_UNDETACHED)
    pthread_attr_setdetachstate(&attr_[i], PTHREAD_CREATE_UNDETACHED);
#elif defined(PTHREAD_CREATE_JOINABLE)
    pthread_attr_setdetachstate(&attr_[i], PTHREAD_CREATE_JOINABLE);
#endif
#ifdef HAVE_PTHREAD_ATTR_GETSTACKSIZE
    size_t defstacksize;
    pthread_attr_getstacksize(&attr_[i], &defstacksize);
#elif HAVE_PTHREAD_ATTR_SETSTACKSIZE
    size_t defstacksize = 1;
#endif
#ifdef HAVE_PTHREAD_ATTR_SETSTACKSIZE
    size_t minstacksize = 2097152;
    if (defstacksize < minstacksize) {
      pthread_attr_setstacksize(&attr_[i], minstacksize);
    }
#endif
  }
}

void PthreadThreadGrp::add_thread(int ithread, Thread* t, int priority)
{
  if (ithread >= nthread_) {
    ExEnv::err0() << indent
                 << "PthreadThreadGrp::add_thread(int, Thread*, int, int): trying to"
                 << "add too many threads" << endl;
  }
  else {
    threads_[ithread] = t;
    //init_priority(ithread, priority);
  }
  
}

#if defined(HAVE_SCHED_GET_PRIORITY_MAX) \
   && defined(HAVE_SCHED_GET_PRIORITY_MIN) \
   && defined(HAVE_PTHREAD_ATTR_SETSCOPE) \
   && defined(HAVE_PTHREAD_ATTR_SETSCHEDPARAM) \
   && defined(HAVE_PTHREAD_ATTR_SETINHERITSCHED) \
   && defined(HAVE_PTHREAD_ATTR_SETSCHEDPOLICY)
#define THREAD_PRIORITY_CAN_BE_SET
#else
#undef THREAD_PRIORITY_CAN_BE_SET
#endif

void PthreadThreadGrp::init_priority(int ithread, int priority)
{
#ifdef THREAD_PRIORITY_CAN_BE_SET
  struct sched_param param, low_param, high_param;
  int rc, selected_sched, set_params;

  set_params=0;
  
  // Check priority settings for various schedulers and select which to use
  selected_sched=-1;

#ifdef SCHED_OTHER
  low_param.sched_priority = sched_get_priority_min(SCHED_OTHER);
  high_param.sched_priority = sched_get_priority_max(SCHED_OTHER);
  if (high_param.sched_priority > low_param.sched_priority) {
    selected_sched = SCHED_OTHER;
    set_params=1;
  }
#endif // SCHED_OTHER
#ifdef SCHED_RR
  if (!set_params) {
    low_param.sched_priority = sched_get_priority_min(SCHED_RR);
    high_param.sched_priority = sched_get_priority_max(SCHED_RR);
    if (high_param.sched_priority > low_param.sched_priority) {
      selected_sched=SCHED_RR; set_params=1;
    }
  }
#endif // SCHED_RR
#ifdef SCHED_FIFO
  if (!set_params) {
    low_param.sched_priority = sched_get_priority_min(SCHED_FIFO);
    high_param.sched_priority = sched_get_priority_max(SCHED_FIFO);
    if (high_param.sched_priority > low_param.sched_priority) {
      selected_sched=SCHED_FIFO; set_params=1;
    }
  }
#endif // SCHED_FIFO

#ifdef PTHREAD_SCOPE_SYSTEM
  pthread_attr_setscope(&attr_[ithread],PTHREAD_SCOPE_SYSTEM);
#endif
  if (set_params) {  
    pthread_attr_setinheritsched(&attr_[ithread],PTHREAD_EXPLICIT_SCHED);
    pthread_attr_setschedpolicy(&attr_[ithread], selected_sched);
    param.sched_priority = ( sched_get_priority_min(selected_sched) + priority );
    pthread_attr_setschedparam(&attr_[ithread],&param);
  }
#endif // THREAD_PRIORITY_CAN_BE_SET
}
int
PthreadThreadGrp::start_threads()
{
  for (int i=1; i < nthread_; i++) {
    if (threads_[i]) {
      int res = pthread_create(&pthreads_[i], &attr_[i],
                               Thread__run_Thread_run,
                               (void*) threads_[i]);
      if (res) {
        ExEnv::errn() << indent << "pthread_create failed" << endl;
        ExEnv::errn() << "error: " << res << ": " << strerror(res) << endl;
        return -1;
      }
    }
  }
  if (threads_[0]) threads_[0]->run();

  return 0;
}

int
PthreadThreadGrp::wait_threads()
{
  for (int i=1; i < nthread_; i++) {
    void *tn;
    if (threads_[i]) {
      int rc = pthread_join(pthreads_[i], (void**)&tn);
      if (rc) {
        ExEnv::errn()
          << "PthreadThreadGrp::wait_threads(): error joining thread"
          << endl;
        ExEnv::errn() << "error: " << rc << ": " << strerror(rc) << endl;
        abort();
      }
    }
  }
    
  return 0;
}

Ref<ThreadLock>
PthreadThreadGrp::new_lock()
{
  return new PthreadThreadLock;
}

ThreadGrp*
PthreadThreadGrp::clone(int nthread)
{
  return new PthreadThreadGrp(*this,nthread);
}

/////////////////////////////////////////////////////////////////////////////

}

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
