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

#ifdef __GNUC__
#pragma implementation
#endif

extern "C" {
#include <userp.h>
#include <lock.h>
}

#include <util/keyval/keyval.h>
#include <util/group/thpuma.h>
#include <util/misc/formio.h>

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

#define CLASSNAME PumaThreadGrp
#define PARENTS public ThreadGrp
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
PumaThreadGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = ThreadGrp::_castdown(cd);
  return do_castdowns(casts,cd);
}

PumaThreadGrp::PumaThreadGrp()
  : ThreadGrp()
{
}

PumaThreadGrp::PumaThreadGrp(const RefKeyVal& keyval)
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

int
PumaThreadGrp::start_threads()
{
  flag_=0;
  if (nthread_ > 1)
    cop((void (*)(void*))threads_[1]->run, &flag, (void*)threads_[1]);

  threads_[0]->run(threads_[0]);
  
  return 0;
}

int
PumaThreadGrp::wait_threads()
{
  if (nthread_ > 1)
    while (!flag);
    
  return 0;
}

RefThreadLock
PumaThreadGrp::new_lock()
{
  return new PumaThreadLock;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
