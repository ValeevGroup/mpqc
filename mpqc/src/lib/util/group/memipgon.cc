//
// memipgon.cc
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

#ifndef _util_group_memipgon_cc
#define _util_group_memipgon_cc

#ifdef __GNUC__
#pragma implementation
#endif

#include <unistd.h>
#include <util/group/memipgon.h>
#include <util/misc/formio.h>

extern "C" {
#include <nx.h>
void msgwait(long);
}

///////////////////////////////////////////////////////////////////////
// The IParagonMemoryGrp class

#define CLASSNAME IParagonMemoryGrp
#define HAVE_KEYVAL_CTOR
#define PARENTS public MIDMemoryGrp
#include <util/class/classi.h>
void *
IParagonMemoryGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  MIDMemoryGrp::_castdown(cd);
  return do_castdowns(casts,cd);
}

IParagonMemoryGrp::IParagonMemoryGrp(const RefMessageGrp& msg):
  MIDMemoryGrp(msg)
{
  use_acknowledgments_ = 0;
  use_active_messages_ = 0;
}

IParagonMemoryGrp::IParagonMemoryGrp(const RefKeyVal& keyval):
  MIDMemoryGrp(keyval)
{
  use_acknowledgments_ = 0;
  use_active_messages_ = 0;
}

IParagonMemoryGrp::~IParagonMemoryGrp()
{
}

long
IParagonMemoryGrp::lockcomm()
{
  return 0;
}

void
IParagonMemoryGrp::unlockcomm(long oldvalue)
{
}

long
IParagonMemoryGrp::send(void* data, int nbytes, int node, int type)
{
  long mid = isend(type, (char*)data, nbytes, node, 0);
  if (debug_) cout << me() << ": IParagonMemoryGrp::send(void*, "
                   << nbytes << ", "
                   << node << ", "
                   << type << "): mid = " << mid
                   << endl;
  return mid;
}

long
IParagonMemoryGrp::recv(void* data, int nbytes, int node, int type)
{
  long mid = irecvx(type, (char*)data, nbytes, node, 0, msginfo);
  if (debug_) cout << me() << ": IParagonMemoryGrp::recv(void*, "
                   << nbytes << ", "
                   << node << ", "
                   << type << "): mid = " << mid
                   << endl;
  return mid;
}

int
IParagonMemoryGrp::probe(long type)
{
  int ret = iprobe((int)type);

  if (ret < 0) {
      cerr << "IParagonMemoryGrp::probe() failed\n";
      return 0;
    }
  
  return ret;
}

long
IParagonMemoryGrp::postrecv(void *data, int nbytes, int type)
{
  cerr << scprintf("IParagonMemoryGrp::postrecv: not implemented\n");
  sleep(1);
  abort();
  return 0;
}

long
IParagonMemoryGrp::wait(long mid1, long mid2)
{
  if (debug_) {
      cout << me() << ": IParagonMemoryGrp::wait("
           << mid1 << ", "
           << mid2 << ")"
           << endl;
    }

  if (mid1 == -1) {
      cerr << me() << ": IParagonMemoryGrp::wait: mid1 == -1" << endl;
      sleep(1);
      abort();
    }
  else if (mid2 == -1) {
      msgwait(mid1);
      if (debug_)
          cout << me() << ": IParagonMemoryGrp::wait(): got " << mid1 << endl;
      return mid1;
    }
  else {
      while(1) {
          if (msgdone(mid1)) {
              if (debug_)
                  cout << me() << ": IParagonMemoryGrp::wait(): got "
                       << mid1 << endl;
              return mid1;
            }
          if (msgdone(mid2)) {
              if (debug_)
                  cout << me() << ": IParagonMemoryGrp::wait(): got "
                       << mid2 << endl;
              return mid2;
            }
        }
    }
}

#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
