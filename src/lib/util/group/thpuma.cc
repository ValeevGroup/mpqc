//
// thpuma.cc
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

extern "C" {
#include <userp.h>
#include <lock.h>
}

#include <util/keyval/keyval.h>
#include <util/group/thpuma.h>
#include <util/misc/formio.h>

using namespace sc;

/////////////////////////////////////////////////////////////////////////////
// PumaThreadLock class

class PumaThreadLock : public ThreadLock {
  private:
    volatile int lock_;
    
  public:
    PumaThreadLock() : lock_(0) {}
    ~PumaThreadLock() {}

    void lock() { set_lock(&lock_); }
    void unlock() { clr_lock(&lock_); }
};

/////////////////////////////////////////////////////////////////////////////
// PumaThreadGrp members

static ClassDesc PumaThreadGrp_cd(
  typeid(PumaThreadGrp),"PumaThreadGrp",1,"public ThreadGrp",
  0, create<PumaThreadGrp>, 0);

PumaThreadGrp::PumaThreadGrp()
  : ThreadGrp()
{
}

PumaThreadGrp::PumaThreadGrp(const Ref<KeyVal>& keyval)
  : ThreadGrp(keyval)
{
  if (nthread_ > 2) {
    delete[] threads_;
    nthread_ = 2;
    threads_ = new Thread*[nthread_];
  }
}

PumaThreadGrp::~PumaThreadGrp()
{
}

static void
run_Thread_run(void *thread)
{
  Thread::run_Thread_run(thread);
}

int
PumaThreadGrp::start_threads()
{
  flag_=0;
  if (nthread_ > 1 && threads_[1])
    cop(run_Thread_run, &flag_, (void*)threads_[1]);

  if (threads_[0]) threads_[0]->run();
  
  return 0;
}

int
PumaThreadGrp::wait_threads()
{
  if (nthread_ > 1)
    while (!flag_);
    
  return 0;
}

Ref<ThreadLock>
PumaThreadGrp::new_lock()
{
  return new PumaThreadLock;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
