//
// memmpl.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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

#ifndef _util_group_memmpl_cc
#define _util_group_memmpl_cc

#ifdef __GNUC__
#pragma implementation
#endif

#include <unistd.h>
#include <util/misc/formio.h>
#include <util/group/memmpl.h>

#include <mpi.h>
#include <mpproto.h>

///////////////////////////////////////////////////////////////////////
// The handler function and its data

static volatile int global_source, global_type, global_mid;
static MPLMemoryGrp *global_mpl_mem = 0;

static void
mpl_memory_handler(int*msgid_arg)
{
  long lmid = *msgid_arg;
  if (!global_mpl_mem) {
      cerr << scprintf("WARNING: Tried to call mpl_memory_handler"
              " without global_mpl_mem\n");
    }
  else {
      global_mpl_mem->handler(&lmid);
    }
}

///////////////////////////////////////////////////////////////////////
// The MPLMemoryGrp class

#define CLASSNAME MPLMemoryGrp
#define HAVE_KEYVAL_CTOR
#define PARENTS public MIDMemoryGrp
#include <util/class/classi.h>
void *
MPLMemoryGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  MIDMemoryGrp::_castdown(cd);
  return do_castdowns(casts,cd);
}

long
MPLMemoryGrp::lockcomm()
{
  int oldvalue;
  mpc_lockrnc(1, &oldvalue);
  if (debug_) cout << ">>>> mpc_lockrnc(1," << oldvalue << ") (lock)" << endl;
  return oldvalue;
}

void
MPLMemoryGrp::unlockcomm(long oldvalue)
{
  int old = oldvalue;
  mpc_lockrnc(old, &old);
  if (debug_) cout << ">>>> mpc_lockrnc(" << oldvalue
                   << "," << old << ") (unlock)" << endl;
}

long
MPLMemoryGrp::send(void* data, int nbytes, int node, int type)
{
  int mid;
  mpc_send(data, nbytes, node, type, &mid);
  if (debug_) cout << ">>>> mpc_send(,"
                   << nbytes << ","
                   << node << ","
                   << type << ","
                   << mid << ")" << endl;
  return mid;
}

long
MPLMemoryGrp::recv(void* data, int nbytes, int node, int type)
{
  int n;
  if (node == -1) n = DONTCARE;
  else n = node;
  int t = type;
  int mid;
  mpc_recv(data, nbytes, &n, &t, &mid);
  if (debug_) cout << ">>>> mpc_recv(,"
                   << nbytes << ","
                   << n << ","
                   << t << ","
                   << mid << ")" << endl;
  return mid;
}

long
MPLMemoryGrp::postrecv(void *data, int nbytes, int type)
{
  global_type = type;
  global_source = DONTCARE; 
  mpc_rcvncall(data, nbytes,
               (int*)&global_source, (int*)&global_type, (int*)&global_mid,
               mpl_memory_handler);
  if (debug_) cout << ">>>> mpc_rcvncall(,"
                   << nbytes << ","
                   << ","
                   << global_type << ","
                   << global_mid << ","
                   << ")" << endl;
  return global_mid;
}

long
MPLMemoryGrp::wait(long mid1, long mid2)
{
  int imid;
  if (mid2 == -1) imid = (int)mid1;
  else imid = DONTCARE;
  size_t count;
  if (debug_)
      cout << scprintf("MPLMemoryGrp: waiting on %d\n", imid);
  if (mpc_wait(&imid,&count)) {
      cerr << scprintf("MPLMemoryGrp: mpc_wait failed\n");
      sleep(1);
      abort();
    }
  if (debug_) cout << ">>>> mpc_wait("
                   << imid << ","
                   << count << ")" << endl;
  return (long)imid;
}

MPLMemoryGrp::MPLMemoryGrp(const RefMessageGrp& msg):
  MIDMemoryGrp(msg)
{
  if (debug_) 
      cout << scprintf("MPLMemoryGrp entered\n");
  if (global_mpl_mem) {
      cerr << scprintf("MPLMemoryGrp: only one allowed at a time\n");
      sleep(1);
      abort();
    }

  global_mpl_mem = this;

  use_acknowledgments_ = 0;
  use_active_messages_ = 1;

  if (debug_) 
      cout << scprintf("MPLMemoryGrp activating\n");
  activate();
  if (debug_) 
      cout << scprintf("MPLMemoryGrp done\n");
}

MPLMemoryGrp::MPLMemoryGrp(const RefKeyVal& keyval):
  MIDMemoryGrp(keyval)
{
  if (debug_) 
      cout << scprintf("MPLMemoryGrp KeyVal entered\n");
  if (global_mpl_mem) {
      cerr << scprintf("MPLMemoryGrp: only one allowed at a time\n");
      sleep(1);
      abort();
    }

  global_mpl_mem = this;

  if (debug_) 
      cout << scprintf("MPLMemoryGrp activating\n");
  activate();
  if (debug_) 
      cout << scprintf("MPLMemoryGrp done\n");
}

MPLMemoryGrp::~MPLMemoryGrp()
{
  if (debug_) 
      cout << scprintf("MPLMemoryGrp: in DTOR\n");
  deactivate();

  int oldlock = lockcomm();
  global_mpl_mem = 0;
  unlockcomm(oldlock);
}

void
MPLMemoryGrp::deactivate()
{
  if (!global_mpl_mem) return;
  MIDMemoryGrp::deactivate();
}

#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
