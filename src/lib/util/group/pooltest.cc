//
// pooltest.cc
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

#include <stdlib.h>
#include <math.h>

#include <util/group/pool.h>
#include <util/misc/formio.h>

using namespace std;
using namespace sc;

class Double {
  private:
    static Pool* pool_;
    static Double* list;
    Double* next;
    Double* prev;
    double* d;
    int size;
    static void zaplist(Double*);
  public:
    Double(size_t size);
    ~Double();
    void zap();
    static void zapall();
    void clear();
    static void pool(Pool*);
};

Double* Double::list = 0;
Pool* Double::pool_ = 0;

void
Double::clear()
{
  if (d) pool_->release(d);
  d = 0;
  size = 0;
}

void
Double::zapall()
{
  zaplist(list);
}

void
Double::zaplist(Double*l)
{
  for (Double* i=l; i; i = i->next) {
      i->zap();
    }
}

Double::Double(size_t s):
  size(s)
{
  if (!pool_) {
      cerr << scprintf("Double::Double: Pool not initialized\n");
      abort();
    }
  d = pool_->allocate_double(size);
  if (!d) {
      //cout << scprintf("\nDouble::Double allocation of size %d failed\n",
      //                 size);
      cout << "F" << endl;
      size = 0;
    }
  next = list;
  prev = 0;
  list = this;
  if (next) next->prev = this;
}

Double::~Double()
{
  clear();
  if (next) next->prev = prev;
  if (prev) prev->next = next;
  else list = next;
}

void
Double::zap()
{
  int i;
  int* x = (int*)d;
  for (i=0; i<size*2; i++) {
      if (x[i] == PoolData::magic) {
          cerr << scprintf("Double::zap: tried to zap a magic number\n");
          abort();
        }
    }
  for (i=0; i<size; i++) d[i] = 0.0;
}

void
Double::pool(Pool*p)
{
  if (pool_ && list) {
      cerr << scprintf("Double::pool: cannot reinitialize pool\n");
      abort();
    }
  pool_ = p;
}

void test1(Pool*);
void test2(Pool*);

int
main(int argc, char* argv[])
{
  const int poolsize = 4000000;
  Pool *pool = new(malloc(poolsize)) Pool(poolsize);

  Double::pool(pool);

  srand48(100);

  cout << "test1:" << endl;
  test1(pool);
  cout << "test2:" << endl;
  test2(pool);

  return 0;
}


void
test1(Pool*pool)
{
  pool->check();

  pool->print();
  cout << endl;

  Double t1(10);

  Double::zapall();

  pool->check();

  pool->print();
  cout << endl;

  Double t2(10000);

  Double::zapall();

  pool->check();

  Double t3(100);

  Double::zapall();

  pool->check();
  
  pool->print();
  cout << endl;

  Double::zapall();

  pool->check();
  
  pool->print();
  cout << endl;

  t2.clear();

  pool->check();
  
  pool->print();
  cout << endl;

  Double t4(100);

  pool->check();

  pool->print();
  cout << endl;

  t1.clear();
  t4.clear();
  t3.clear();

  pool->check();

  pool->print();
  cout << endl;

}

void
test2(Pool*pool)
{
  int i, ii;

  const int npass = 10;
  const int nd = 512;
  Double* d[nd];

  for (i=0; i<nd; i++) d[i] = 0;
       
  for (ii=0; ii<npass; ii++) {
      for (i=0; i<nd;) {
          if (mrand48() > 0) {
              // allocate data
              size_t size = lrand48() & 0x03ff;
              d[i] = new Double(size);
              cerr.flush();
              cout << "a" << flush;
              i++;
            }
          else {
              // deallocate data
              int loc = (int) (drand48()*i);
              if (loc >= nd) loc = nd-1;
              if (loc < 0) loc = 0;
              if (d[loc]) {
                  cerr.flush();
                  cout << "d" << flush;
                  delete d[loc];
                  d[loc] = 0;
                }
            }
          //pool->print();
          //pool->check();
          //Double::zapall();
        }
      for (i=0; i<nd; i++) {
          if (d[i]) {
              cerr.flush();
              cout << "d" << flush;
              delete d[i];
              d[i] = 0;
            }
        }
      cout << endl;
      pool->print();
      //pool->check();
      //Double::zapall();
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
