//
// rnglock.h
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

#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_rnglock_h
#define _util_group_rnglock_h

#include <iostream>

#include <util/misc/exenv.h>

namespace sc {

class Pool;

class RangeLockItem {
  public:
    RangeLockItem *prev;
    RangeLockItem *next;
    int start;
    int fence;
    int value;
    RangeLockItem(RangeLockItem *p, RangeLockItem *n, int s, int f, int v):
      prev(p), next(n), start(s), fence(f), value(v) {}
    ~RangeLockItem() {};

    static void *operator new(size_t, Pool *);
    static void operator delete(void *, Pool *);
};

class RangeLockValOp;
class RangeLock {
  private:
    RangeLockItem *root_;
    Pool *pool_;

    void split_ranges(int start, int fence);
    void do_valop(RangeLockValOp&, int start, int fence);
  public:
    RangeLock(Pool *pool = 0);
    ~RangeLock();

    void increment(int start, int fence);
    void decrement(int start, int fence);
    void set(int start, int fence, int value);
    void sum(int start, int fence, int delta);

    // check for anything within a range to be equal to a value
    int checkeq(int start, int fence, int value);
    // check for anything within a range to be greater than a value
    int checkgr(int start, int fence, int value);

    void check();
    void print(std::ostream &o = ExEnv::out0()) const;

    int lockvalue(int i);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
