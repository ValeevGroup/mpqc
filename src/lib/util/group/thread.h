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

#include <util/class/class.h>

namespace sc {

/** The ThreadLock abstract class provides mutex locks to be used in
    conjunction with ThreadGrp's.  ThreadLock objects should be
    locked and unlocked with ThreadLockHolder objects to provide
    exception safety.
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


/** Acquire a lock on creation and release it on destruction.
    This should be used to lock and unlock ThreadLock objects
    to provide exception safety.
 */
class ThreadLockHolder {
    Ref<ThreadLock> lock_;
    bool locked_;
  public:
    /// Acquires the lock.
    ThreadLockHolder(const Ref<ThreadLock> &l): lock_(l) {
      lock_->lock();
      locked_ = true;
    }
    /// Release the lock before the DTOR is called, if it is still held.
    void unlock() { if (locked_) { lock_->unlock(); locked_ = false; } }
    /// Acquire the lock once more.
    void lock() { if (!locked_) { lock_->lock(); locked_ = true; } }
    /// Releases the lock if it is still held.
    ~ThreadLockHolder() { unlock(); }
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
    virtual void add_thread(int threadnum, Thread* thread);
    /** Like add_thread(threadnum, thread), but assign a priority that the
        thread is to use.  The member is primarily for experimentation, the
        priority argument is currently not well defined and ignored.  */
    virtual void add_thread(int threadnum, Thread* thread, int priority);
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

    /// Sets the default ThreadGrp. This will be returned
    /// by future calls to get_default_threadgrp.
    static void set_default_threadgrp(const Ref<ThreadGrp>&);
    /// Returns the default ThreadGrp. If set_default_threadgrp
    /// has not already been called to set a default, this will
    /// construct a new ThreadGrp.
    static ThreadGrp * get_default_threadgrp();
    /// Create a ThreadGrp. First, this will determine if the -threadgrp
    /// option has been given with \p argc and \p argv. If so, then
    /// the argument should be ParsedKeyVal input for a ThreadGrp. If the
    /// argument is not found, then the THREADGRP environment variable
    /// is examined. If found, then its value should be a ParsedKeyVal
    /// input for a ThreadGrp. Otherwise null is returned.
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
    
    Ref<ThreadLock> new_lock();

    ThreadGrp* clone(int nthread = -1);
};

}

extern "C" {
    // a C linkage interface to run_Thread_run
    void *Thread__run_Thread_run(void*thread);
}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
