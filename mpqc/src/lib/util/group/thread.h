//
// thread.h
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

#ifndef _util_group_thread_h
#define _util_group_thread_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/class/class.h>

//. The \clsnm{ThreadLock} abstract class provides mutex locks to be used in
//. conjunction with \clsnm{ThreadGrp}'s
class ThreadLock : public VRefCount {
  public:
    ThreadLock();
    virtual ~ThreadLock();

    virtual void lock() =0;
    virtual void unlock() =0;
};
REF_dec(ThreadLock);

//. The \clsnm{Thread} abstract class defines an interface which must be
//. implemented by classes wishing to be run as threads.
class Thread {
  public:
    Thread();
    virtual ~Thread();

    static void *run_Thread_run(void*thread);
    virtual void run() =0;
};
    
DescribedClass_REF_fwddec(ThreadGrp);

//. The \clsnm{ThreadGrp} abstract class privides a means to manage separate
//. threads of control.
class ThreadGrp: public DescribedClass {
#define CLASSNAME ThreadGrp
#include <util/class/classda.h>
  protected:
    Thread** threads_;
    int nthread_;

  public:
    ThreadGrp();
    ThreadGrp(const RefKeyVal&);
    virtual ~ThreadGrp();

    void add_thread(int, Thread*);
    int nthread() const { return nthread_; }

    virtual int start_threads() =0;
    virtual int wait_threads() =0;
    virtual RefThreadLock new_lock() =0;

    static void set_default_threadgrp(const RefThreadGrp&);
    static ThreadGrp * get_default_threadgrp();
    static ThreadGrp * initial_threadgrp(int &argc, char ** argv);
};
DescribedClass_REF_dec(ThreadGrp);

//. The \clsnm{ProcThreadGrp} class privides a concrete thread group
//. appropriate for an environment where there is only one thread.
class ProcThreadGrp: public ThreadGrp {
#define CLASSNAME ProcThreadGrp
#define HAVE_KEYVAL_CTOR
#include <util/class/classd.h>
  public:
    ProcThreadGrp();
    ProcThreadGrp(const RefKeyVal&);
    ~ProcThreadGrp();

    int start_threads();
    int wait_threads();
    RefThreadLock new_lock();
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
