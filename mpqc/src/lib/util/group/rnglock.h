
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_rnglock_h
#define _util_group_rnglock_h

#include <iostream.h>

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

    void *operator new(size_t, Pool *);
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
    void print(ostream &o = cout);

    int lockvalue(int i);
};

#endif
