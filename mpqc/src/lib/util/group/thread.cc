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

#ifdef __GNUC__
#pragma implementation
#endif

#include <util/keyval/keyval.h>
#include <util/group/thread.h>
#include <util/misc/formio.h>
#include <util/misc/exenv.h>

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

#define CLASSNAME ThreadGrp
#define PARENTS public DescribedClass
#include <util/class/classia.h>
void *
ThreadGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

ThreadGrp::ThreadGrp() : threads_(0),  nthread_(1)
{
  threads_ = new Thread*[nthread_];
}

ThreadGrp::ThreadGrp(const ThreadGrp &tg, int nthread)
{
  if (nthread == -1) nthread_ = tg.nthread_;
  else nthread_ = nthread;
  threads_ = new Thread*[nthread_];
}

ThreadGrp::ThreadGrp(const RefKeyVal& keyval)
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
    cerr << node0 << indent
         << "ThreadGrp::add_thread: trying to add too many threads"
         << endl;
  } else {
    threads_[i] = t;
  }
}

static RefThreadGrp default_threadgrp;

void
ThreadGrp::set_default_threadgrp(const RefThreadGrp& grp)
{
  default_threadgrp = grp;
}

ThreadGrp*
ThreadGrp::get_default_threadgrp()
{
  if (default_threadgrp.null())
    default_threadgrp = new ProcThreadGrp;

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
        i++;
        if (i >= argc) {
          cerr << "-threadgrp must be following by an argument"
               << endl;
          abort();
        }
        keyval_string = argv[i];
        // permute the messagegrp arguments to the end of argv
        char *tmp = argv[argc-2];
        argv[argc-2] = argv[i-1];
        argv[i-1] = tmp;
        tmp = argv[argc-1];
        argv[argc-1] = argv[i];
        argv[i] = tmp;
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
    RefParsedKeyVal strkv = new ParsedKeyVal();
    strkv->parse_string(keyval_string);
    RefDescribedClass dc = strkv->describedclassvalue();
    grp = ThreadGrp::castdown(dc.pointer());
    if (dc.null()) {
      cerr << "initial_threadgrp: couldn't find a ThreadGrp in "
           << keyval_string << endl;
      abort();
    } else if (!grp) {
      cerr << "initial_threadgrp: wanted ThreadGrp but got "
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
  cout << "ThreadGrp::clone not supported for " << class_name() << endl;
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

#define CLASSNAME ProcThreadGrp
#define PARENTS public ThreadGrp
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
ProcThreadGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = ThreadGrp::_castdown(cd);
  return do_castdowns(casts,cd);
}

ProcThreadGrp::ProcThreadGrp()
  : ThreadGrp()
{
}

ProcThreadGrp::ProcThreadGrp(const RefKeyVal& keyval)
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

RefThreadLock
ProcThreadGrp::new_lock()
{
  return new ProcThreadLock;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
