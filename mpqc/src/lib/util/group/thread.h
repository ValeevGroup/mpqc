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

/** The ThreadLock abstract class provides mutex locks to be used in
    conjunction with ThreadGrp's.
*/
class ThreadLock : public RefCount {
  public:
    ThreadLock();
    virtual ~ThreadLock();

    /// Obtain the lock.
    virtual void lock() =0;
    /// Release the lock.
    virtual void unlock() =0;
};


/** The Thread abstract class defines an interface which must be
    implemented by classes wishing to be run as threads. */
class Thread {
  public:
    Thread();
    virtual ~Thread();

    static void *run_Thread_run(void*thread);

    /// This is called with the Thread is run from a ThreadGrp.
    virtual void run() =0;
};

extern "C" {
    // a C linkage interface to run_Thread_run
    void *Thread__run_Thread_run(void*thread);
}
    
/** The ThreadGrp abstract class provides a means to manage separate
    threads of control. */
class ThreadGrp: public DescribedClass {
  protected:
    Thread** threads_;
    int nthread_;

  public:
    ThreadGrp();
    ThreadGrp(const Ref<KeyVal>&);
    ThreadGrp(const ThreadGrp&, int nthread = -1);
    virtual ~ThreadGrp();

    /** Assigns a Thread object to each thread.  If 0 is assigned to
        a thread, then that thread will be skipped. */
    virtual void add_thread(int, Thread*);
    virtual void add_thread(int, Thread*, int);
    /// The number of threads that will be run by start_thread.
    int nthread() const { return nthread_; }

    void delete_threads();

    /** Starts the threads running.  Thread 0 will be run by the
        thread that calls start_threads. */
    virtual int start_threads() =0;
    /** Wait for all the threads to complete.  This must be called
        before start_threads is called again or the object is destroyed. */
    virtual int wait_threads() =0;
    /// Return a local object.
    virtual Ref<ThreadLock> new_lock() =0;

    /** Create a ThreadGrp like the current one.  If nthread is given, the
        new ThreadGrp will attempt to support that number of threads, but
        the actual number supported may be less.  If nthread is -1, the
        number of threads in the current group will be used. */
    virtual ThreadGrp* clone(int nthread = -1);

    static void set_default_threadgrp(const Ref<ThreadGrp>&);
    static ThreadGrp * get_default_threadgrp();
    static ThreadGrp * initial_threadgrp(int &argc, char ** argv);
};


/** The ProcThreadGrp class privides a concrete thread group
    appropriate for an environment where there is only one thread.
*/
class ProcThreadGrp: public ThreadGrp {
  public:
    ProcThreadGrp();
    ProcThreadGrp(const Ref<KeyVal>&);
    ~ProcThreadGrp();

    int start_threads();
    int wait_threads();
//    void add_thread(int, Thread*);
    void add_thread(int i, Thread*t, int sched, int priority) { ThreadGrp::add_thread(i,t); }
    
    Ref<ThreadLock> new_lock();

    ThreadGrp* clone(int nthread = -1);
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
