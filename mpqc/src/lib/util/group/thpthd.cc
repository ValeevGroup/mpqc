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

#ifdef __GNUC__
#pragma implementation
#endif

#ifdef HAVE_CONFIG_H
#include <scconfig.h>
#endif

#if HAVE_PTHREAD_H
#include <pthread.h>
#endif

#include <util/keyval/keyval.h>
#include <util/group/thpthd.h>
#include <util/misc/formio.h>

using namespace std;

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

#define CLASSNAME PthreadThreadGrp
#define PARENTS public ThreadGrp
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
PthreadThreadGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = ThreadGrp::_castdown(cd);
  return do_castdowns(casts,cd);
}

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

PthreadThreadGrp::PthreadThreadGrp(const RefKeyVal& keyval)
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
  }
  delete attr_;
}

void
PthreadThreadGrp::init_attr()
{
  attr_ = new pthread_attr_t;
  pthread_attr_init(attr_);
#if defined(PTHREAD_CREATE_UNDETACHED)
  pthread_attr_setdetachstate(attr_, PTHREAD_CREATE_UNDETACHED);
#elif defined(PTHREAD_CREATE_JOINABLE)
  pthread_attr_setdetachstate(attr_, PTHREAD_CREATE_JOINABLE);
#endif
#ifdef HAVE_PTHREAD_ATTR_GETSTACKSIZE
  size_t defstacksize;
  pthread_attr_getstacksize(attr_, &defstacksize);
#elif HAVE_PTHREAD_ATTR_SETSTACKSIZE
  size_t defstacksize = 1;
#endif
#ifdef HAVE_PTHREAD_ATTR_SETSTACKSIZE
  size_t minstacksize = 2097152;
  if (defstacksize < minstacksize) {
    pthread_attr_setstacksize(attr_, minstacksize);
  }
#endif
}

int
PthreadThreadGrp::start_threads()
{
  for (int i=1; i < nthread_; i++) {
    if (threads_[i]) {
      int res = pthread_create(&pthreads_[i], attr_,
                               Thread::run_Thread_run,
                               (void*) threads_[i]);
      if (res) {
        ExEnv::err() << indent << "thread death " << res << endl;
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
    int tn;
    if (threads_[i]) {
      if (pthread_join(pthreads_[i], (void**)&tn)) {
        ExEnv::out()
          << "PthreadThreadGrp::wait_threads(): error joining thread"
          << endl;
        abort();
      }
    }
  }
    
  return 0;
}

RefThreadLock
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

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
