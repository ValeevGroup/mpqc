//
// thrtest.cc
// based on: messtest.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
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
#include <util/misc/formio.h>
#include <util/misc/bug.h>
#include <util/group/thread.h>

using namespace std;
using namespace sc;

// Force linkages:
//#ifndef __PIC__
#ifdef PUMAGON
#   include <util/group/thpuma.h>
    static ForceLink<PumaThreadGrp> fl0;
#endif
# ifdef HAVE_PTHREAD
#   include <util/group/thpthd.h>
    static ForceLink<PthreadThreadGrp> fl2;
# endif
//#endif

class TestThread: public Thread {
  private:
    Ref<ThreadLock> lock;
  public:
    static int count;
    TestThread(const Ref<ThreadLock> &l): lock(l) {}
    void run();
    int n() const { return 1000000; }
};

int TestThread::count = 0;

void
TestThread::run()
{
  for (int i=0; i<n(); i++) {
      lock->lock();
      count++;
      lock->unlock();
    }
}

int
main(int argc, char**argv)
{
  int i;

  Ref<ThreadGrp> grp = ThreadGrp::initial_threadgrp(argc, argv);

  Ref<Debugger> debugger;

  if (grp == 0) {
      const char* input = SRCDIR "/thrtest.in";
      const char* keyword = "thread";

      if (argc >= 2) input = argv[1];
      if (argc >= 3) keyword = argv[2];

      Ref<KeyVal> keyval = new ParsedKeyVal(input);

      grp << keyval->describedclassvalue(keyword);

      debugger << keyval->describedclassvalue(":debug");

      if (grp == 0) {
          cerr << scprintf("Couldn't initialize ThreadGrp\n");
          abort();
        }
    }

  if (debugger) {
      debugger->set_exec(argv[0]);
      debugger->set_prefix(0);
    }

  Debugger::set_default_debugger(debugger);

  TestThread **thr = new TestThread*[grp->nthread()];
  Ref<ThreadLock> lock = grp->new_lock();
  for (i=0; i<grp->nthread(); i++) {
      thr[i] = new TestThread(lock);
      grp->add_thread(i,thr[i]);
    }

  grp->start_threads();
  grp->wait_threads();

  int right_count = 0;
  for (i=0; i<grp->nthread(); i++) {
      right_count += thr[i]->n();
      delete thr[i];
    }
  delete[] thr;

  cout << "count = " << TestThread::count << " (expected " << right_count << ")" << endl;

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
