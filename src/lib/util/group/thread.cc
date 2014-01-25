//
// thread.cc
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

#include <util/keyval/keyval.h>
#include <util/group/thread.h>
#include <util/misc/formio.h>
#include <util/misc/exenv.h>

// debug includes
#include <sys/types.h>
#include <unistd.h>

#include <mpqc_config.h>
#ifdef HAVE_PTHREAD
#  include <util/group/thpthd.h>
#endif

using namespace std;
using namespace sc;

// This has C linkage.
void *
Thread__run_Thread_run(void* vth)
{
  return Thread::run_Thread_run(vth);
}

namespace sc {

/////////////////////////////////////////////////////////////////////////////

ThreadLock::ThreadLock()
{
}

ThreadLock::~ThreadLock()
{
}

/////////////////////////////////////////////////////////////////////////////

Thread::Thread()
{
}

Thread::~Thread()
{
}

void *
Thread::run_Thread_run(void* vth)
{
  if (vth) ((Thread*)vth)->run();
  return 0;
}

/////////////////////////////////////////////////////////////////////////////
// ThreadGrp members

static ClassDesc ThreadGrp_cd(typeid(ThreadGrp),"ThreadGrp",1,
                              "public DescribedClass");

ThreadGrp::ThreadGrp() : threads_(0),  nthread_(1)
{
  threads_ = new Thread*[nthread_];
  for (int i = 0; i<nthread_; i++) threads_[i] = 0;
}

ThreadGrp::ThreadGrp(const ThreadGrp &tg, int nthread)
{
  if (nthread == -1) nthread_ = tg.nthread_;
  else nthread_ = nthread;
  threads_ = new Thread*[nthread_];
  for (int i = 0; i<nthread_; i++) threads_[i] = 0;
}

ThreadGrp::ThreadGrp(const Ref<KeyVal>& keyval)
{
  int defaultnum = ExEnv::nproc();
  if (defaultnum == 0) defaultnum = 1;
  KeyValValueint num(defaultnum);
  nthread_ = keyval->intvalue("num_threads",num);

  threads_ = new Thread*[nthread_];
  for (int i=0; i<nthread_; i++) threads_[i] = 0;
}

ThreadGrp::~ThreadGrp()
{
  if (nthread_) {
    delete[] threads_;
    nthread_=0;
    threads_=0;
  }
}

void
ThreadGrp::delete_threads()
{
  for (int i=0; i<nthread_; i++) {
    delete threads_[i];
    threads_[i] = 0;
  }
}

void
ThreadGrp::add_thread(int i, Thread*t)
{
  if (i >= nthread_) {
    ExEnv::err0() << indent
         << "ThreadGrp::add_thread: trying to add too many threads"
         << endl;
  } else {
    threads_[i] = t;
  }
}

void
ThreadGrp::add_thread(int i, Thread*t, int priority)
{
  add_thread(i,t);
}

static Ref<ThreadGrp> default_threadgrp;

void
ThreadGrp::set_default_threadgrp(const Ref<ThreadGrp>& grp)
{
  default_threadgrp = grp;
}

ThreadGrp*
ThreadGrp::get_default_threadgrp()
{
  if (default_threadgrp.null()) {
#ifdef HAVE_PTHREAD
    default_threadgrp = new PthreadThreadGrp;
#else
    default_threadgrp = new ProcThreadGrp;
#endif
  }

  return default_threadgrp;
}

ThreadGrp*
ThreadGrp::initial_threadgrp(int& argc, char ** argv)
{
  ThreadGrp *grp = 0;
  char * keyval_string = 0;
  
  // see if a thread group is given on the command line
  if (argc && argv) {
    for (int i=0; i < argc; i++) {
      if (argv[i] && !strcmp(argv[i], "-threadgrp")) {
        char *threadgrp_string = argv[i];
        i++;
        if (i >= argc) {
          ExEnv::errn() << "-threadgrp must be following by an argument"
               << endl;
          abort();
        }
        keyval_string = argv[i];
        // move the threadgrp arguments to the end of argv
        int j;
        for (j=i+1; j<argc; j++) {
          argv[j-2] = argv[j];
        }
        argv[j++] = threadgrp_string;
        argv[j++] = keyval_string;
        // decrement argc to hide the last two arguments
        argc -= 2;
        break;
      }
    }
  }

  if (!keyval_string) {
    // find out if the environment gives the containing thread group
    keyval_string = getenv("THREADGRP");
    if (keyval_string) {
      if (!strncmp("THREADGRP=", keyval_string, 11)) {
        keyval_string = strchr(keyval_string, '=');
      }
      if (*keyval_string == '=') keyval_string++;
    }
  }

  // if keyval input for a thread group was found, then
  // create it.
  if (keyval_string) {
    if (keyval_string[0] == '\0') return 0;
    Ref<ParsedKeyVal> strkv = new ParsedKeyVal();
    strkv->parse_string(keyval_string);
    Ref<DescribedClass> dc = strkv->describedclassvalue();
    grp = dynamic_cast<ThreadGrp*>(dc.pointer());
    if (dc.null()) {
      ExEnv::errn() << "initial_threadgrp: couldn't find a ThreadGrp in "
           << keyval_string << endl;
      abort();
    } else if (!grp) {
      ExEnv::errn() << "initial_threadgrp: wanted ThreadGrp but got "
           << dc->class_name() << endl;
      abort();
    }
    // prevent an accidental delete
    grp->reference();
    strkv = 0;
    dc = 0;
    // accidental delete not a problem anymore since all smart pointers
    // to grp are dead
    grp->dereference();
    return grp;
  }

  return 0;
}

ThreadGrp*
ThreadGrp::clone(int nthread)
{
  ExEnv::errn() << "ThreadGrp::clone not supported for " << class_name()
               << endl;
  abort();
  return 0;
}

/////////////////////////////////////////////////////////////////////////////
// ProcThreadLock class

class ProcThreadLock : public ThreadLock {
  public:
    ProcThreadLock() {}
    ~ProcThreadLock() {}

    void lock() {}
    void unlock() {}
};

/////////////////////////////////////////////////////////////////////////////
// ProcThreadGrp members

static ClassDesc ProcThreadGrp_cd(
  typeid(ProcThreadGrp),"ProcThreadGrp",1,"public ThreadGrp",
  0, create<ProcThreadGrp>, 0);

ProcThreadGrp::ProcThreadGrp()
  : ThreadGrp()
{
}

ProcThreadGrp::ProcThreadGrp(const Ref<KeyVal>& keyval)
  : ThreadGrp(keyval)
{
  if (nthread_ > 1) {
    delete[] threads_;
    nthread_ = 1;
    threads_ = new Thread*[nthread_];
  }
}

ProcThreadGrp::~ProcThreadGrp()
{
}

int
ProcThreadGrp::start_threads()
{
  if (threads_[0]) threads_[0]->run();
  return 0;
}

int
ProcThreadGrp::wait_threads()
{
  return 0;
}

Ref<ThreadLock>
ProcThreadGrp::new_lock()
{
  return new ProcThreadLock;
}

ThreadGrp*
ProcThreadGrp::clone(int nthread)
{
  return new ProcThreadGrp;
}

/////////////////////////////////////////////////////////////////////////////

}

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
