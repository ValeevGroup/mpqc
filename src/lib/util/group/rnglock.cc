//
// rnglock.cc
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

#include <new>

#include <util/misc/formio.h>
#include <util/group/rnglock.h>
#include <util/group/pool.h>

using namespace std;
using namespace sc;

#define CHECK 0

namespace sc {

/////////////////////////////////////////////////////////////////////
// RangeLockItem members

void *
RangeLockItem::operator new(size_t size, Pool * pool)
{
  if (pool) return pool->allocate(size);
  else return ::operator new(size);
}

void
RangeLockItem::operator delete(void* r, Pool *pool)
{
 if (pool) pool->release(r);
 else ::operator delete(r);
}

/////////////////////////////////////////////////////////////////////
// Utility classes

class RangeLockValOp {
  public:
    virtual ~RangeLockValOp() {};
    virtual void op(int&) = 0;
};

class RangeLockValSum: public RangeLockValOp {
  private:
    int delta;
  public:
    RangeLockValSum(int d): delta(d) {}
    void op(int &i) { i += delta; }
};

class RangeLockValSet: public RangeLockValOp {
  private:
    int val;
  public:
    RangeLockValSet(int v): val(v) {}
    void op(int &i) { i = val; }
};

/////////////////////////////////////////////////////////////////////
// Members of RangeLock

RangeLock::RangeLock(Pool *pool)
{
  pool_ = pool;
  root_ = 0;
}

RangeLock::~RangeLock()
{
  for (RangeLockItem *i = root_; i;) {
      RangeLockItem *next = i->next;
      RangeLockItem::operator delete(i, pool_);
      i = next;
    }
}

int
RangeLock::lockvalue(int loc)
{
  for (RangeLockItem *i = root_; i; i = i->next) {
      if (loc >= i->start && loc < i->fence) {
          return i->value;
        }
      if (loc < i->fence) break;
    }
  return 0;
}

int
RangeLock::checkeq(int start, int fence, int value)
{
  split_ranges(start, fence);
  for (RangeLockItem *i = root_; i; i = i->next) {
      if (((start >= i->start) && (start < i->fence))
          || ((fence > i->start) && (fence <= i->fence))
          || ((start < i->start) && (fence > i->fence))) {
          if (value != i->value) return 0;
        }
      if (fence < i->fence) break;
    }
  return 1;
}

int
RangeLock::checkgr(int start, int fence, int value)
{
  split_ranges(start, fence);
  for (RangeLockItem *i = root_; i; i = i->next) {
      if (((start >= i->start) && (start < i->fence))
          || ((fence > i->start) && (fence <= i->fence))
          || ((start < i->start) && (fence > i->fence))) {
          if (i->value <= value) return 0;
        }
      if (fence < i->fence) break;
    }
  return 1;
}

void
RangeLock::check()
{
  for (RangeLockItem* i = root_; i; i = i->next) {
      int bad = 0;
      if (i->next && i->next->prev != i) {
          ExEnv::errn() << scprintf("i->next->prev bad\n");
          bad = 1;
        }
      if (i->prev && i->prev->next != i) {
          ExEnv::errn() << scprintf("i->prev->next bad\n");
          bad = 1;
        }
      if (i->start >= i->fence) {
          ExEnv::errn() << scprintf("start >= fence\n");
          bad = 1;
        }
      if (i->next && i->fence > i->next->start) {
          ExEnv::errn() << scprintf("fence > next start\n");
          bad = 1;
        }
#if VERBOSE
      ExEnv::outn()
          << scprintf("i = 0x%08x, n = 0x%08x, p = 0x%08x, [%3d, %3d), %5d\n",
                      i, i->next, i->prev, i->start, i->fence, i->value);
#endif
      if (bad) abort();
    }
}

void
RangeLock::split_ranges(int start, int fence)
{
  if (root_ == 0) {
      // no blocks are allocted yet, initialize one and return
      root_ = new(pool_) RangeLockItem(0, 0, start, fence, 0);
      return;
    }

  RangeLockItem *i;
  for (i = root_; i; i=i->next) {
      if (start > i->start && start < i->fence) {
          // start is in the middle of this block, split it
          RangeLockItem *t = new(pool_) RangeLockItem(i, i->next,
                                                      start, i->fence,
                                                      i->value);
          i->fence = start;
          i->next = t;
          if (t->next) t->next->prev = t;
          i = t;
          break;
        }
      else if (start < i->start && fence <= i->start) {
          // start and end are before this block, insert it and return
          RangeLockItem *t = new(pool_) RangeLockItem(i->prev, i,
                                                      start, fence, 0);
          i->prev = t;
          if (t->prev) t->prev->next = t;
          else root_ = t;
          return;
        }
      else if (start < i->start) {
          // start is before this block, fill in the gap
          RangeLockItem *t = new(pool_) RangeLockItem(i->prev, i,
                                                      start, i->start, 0);
          i->prev = t;
          if (t->prev) t->prev->next = t;
          else root_ = t;
          break;
        }
      else if (start == i->start) {
          // start coincides with this block's start
          break;
        }
      else if (i->next == 0) {
          // start is after the last block, make the block and return
          RangeLockItem *t = new(pool_) RangeLockItem(i, 0, start, fence, 0);
          i->next = t;
          return;
        }
      // otherwise start is after this block, continue
    }

  for (; i; i=i->next) {
      if (fence > i->start && i->prev && i->prev->fence < i->start) {
          // fence is after this block and there is a gap before
          RangeLockItem *t = new(pool_) RangeLockItem(i->prev, i,
                                                      i->prev->fence, i->start,
                                                      0);
          i->prev = t;
          if (t->prev) t->prev->next = t;
          else root_ = t;
          i = t;
        }
      else if (fence > i->start && fence < i->fence) {
          // fence is in the middle of this block, split it and return
          RangeLockItem *t = new(pool_) RangeLockItem(i->prev, i,
                                                      i->start, fence,
                                                      i->value);
          i->start = fence;
          i->prev = t;
          if (t->prev) t->prev->next = t;
          else root_ = t;
          return;
        }
      else if (fence > i->fence && !i->next) {
          // fence is after this block and no blocks follow
          RangeLockItem *t = new(pool_) RangeLockItem(i, 0,
                                                      i->fence, fence, 0);
          i->next = t;
          return;
        }
      else if (fence <= i->start) {
          // fence is before this block, fill in the gap
          RangeLockItem *p = i->prev;
          RangeLockItem *t = new(pool_) RangeLockItem(p, i,
                                                      p->fence, fence, 0);
          p->next = t;
          i->prev = t;
          return;
        }
      else if (fence == i->fence) {
          // fence coincides with this block's fence
          return;
        }
      // otherwise fence is after this block, continue
    }

  ExEnv::errn() << scprintf("RangeLock::split_ranges(): got to end\n");
  abort();
}

void
RangeLock::do_valop(RangeLockValOp& op, int start, int fence)
{
  if (start == fence) return;

  split_ranges(start, fence);

#if CHECK
  check();
#endif

  for (RangeLockItem *i = root_; i; i = i->next) {
      if ((start == i->start)
          || (fence == i->fence)
          || ((start < i->start) && (fence > i->fence))) {
          op.op(i->value);
        }
      if (fence < i->fence) break;
    }
}

void
RangeLock::sum(int start, int fence, int delta)
{
  RangeLockValSum sum(delta);
  do_valop(sum, start, fence);
}

void
RangeLock::increment(int start, int fence)
{
  RangeLockValSum sum(1);
  do_valop(sum, start, fence);
}

void
RangeLock::decrement(int start, int fence)
{
  RangeLockValSum sum(-1);
  do_valop(sum, start, fence);
}

void
RangeLock::set(int start, int fence, int value)
{
  RangeLockValSet set(value);
  do_valop(set, start, fence);
}

void
RangeLock::print(ostream &o) const
{
  for (RangeLockItem *i = root_; i; i = i->next) {
      o << scprintf("  RangeLockItem: [%5d, %5d): %4d\n",
              i->start, i->fence, i->value);
    }
}

/////////////////////////////////////////////////////////////////////////////

}

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
