//
// rnglocktest.cc
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

#define TEST 1
#define VERBOSE 0

#include <math.h>
#include <stdlib.h>
#include <util/misc/formio.h>
#include <util/group/rnglock.h>

using namespace std;
using namespace sc;

main()
{
  const int size = 10000;
  int bin[size];
  int i;

  RangeLock lock;

  for (i=0; i<size; i++) bin[i] = 0;

  for (i=0; i<10000; i++) {
      int start = random()%size;
      int length = random()%30;
      if (length == 0) length++;
      int fence = start + length;
      if (fence > size) fence = size;
      int val = random()%2 ? -1:1;
      val = 1;
#if VERBOSE
      ExEnv::out0() << scprintf("adding block %d: start = %d, fence = %d, val = %d\n",
                       i, start, fence, val);
#endif
      lock.sum(start, fence, val);
      for (int j=start; j<fence; j++) {
          bin[j] += val;
        }
    }

  ExEnv::out0() << "printing nonzero deltas" << endl;

  for (i=0; i<size; i++) {
      int delta = bin[i] - lock.lockvalue(i);
      if (delta) ExEnv::out0() << scprintf(" %5d: %8d (%8d %8d)\n", i, delta,
                                  bin[i], lock.lockvalue(i));
    }

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
