//
// reftest.cc
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

#include <util/ref/ref.h>
#include <util/ref/reftestx.h>

#if HAVE_PTHREAD==1 && REF_USE_LOCKS==1
#include <pthread.h>
#endif

using namespace std;

using sc::Y;
using sc::X;
using sc::Ref;

void
test1()
{
  Y *y1p = new Y;
  Y *y2p = new Y;
  X *x1p = y1p;
  X *x2p = y2p;
  void *vy1p = (void*) y1p;
  void *vy2p = (void*) y2p;
  void *vx1p = (void*) x1p;
  void *vx2p = (void*) x2p;
  Ref<Y> y1 = y1p;
  Ref<Y> y2 = y2p;
  cout << "X::nx after Ref<Y> assignment: " << X::nx << "(2)" << endl;
  cout << "y1->nreference() after Ref<Y> assignment: "
       << y1->nreference() << "(1)" << endl;
  Ref<X> x1(y1);
  Ref<X> x2(y2);
  cout << "X::nx after Ref<X> assignment: " << X::nx << "(2)" << endl;
  cout << "y1->nreference() after Ref<X> assignment: "
       << y1->nreference() << "(2)" << endl;
  x1 = y1;
  x2 = y2;
  Ref<Y> yb1, yb2;
  yb1 << x1;
  yb2 << x2;

  cout << "x1 = " << (void*)x1.pointer() << "(" << vx1p << ")" << endl;
  cout << "x2 = " << (void*)x2.pointer() << "(" << vx2p << ")" << endl;
  cout << "y1 = " << (void*)y1.pointer() << "(" << vy1p << ")" << endl;
  cout << "y2 = " << (void*)y2.pointer() << "(" << vy2p << ")" << endl;
  cout << "yb1 = " << (void*)yb1.pointer()
       << "(" << (void*)y1.pointer() << ")" << endl;
  cout << "yb2 = " << (void*)yb2.pointer()
       << "(" << (void*)y2.pointer() << ")" << endl;
}

void
test2()
{
  Ref<Y> y1 = new Y;
  Ref<Y> y2 = new Y;
  cout << "X::nx after Ref<Y> assignment: " << X::nx << "(2)" << endl;
  cout << "y1->nreference() after Ref<Y> assignment: "
       << y1->nreference() << "(1)" << endl;
  Ref<X> x1(y1);
  Ref<X> x2(y2);
  cout << "X::nx after Ref<X> assignment: " << X::nx << "(2)" << endl;
  cout << "y1->nreference() after Ref<X> assignment: "
       << y1->nreference() << "(2)" << endl;
  x1 = y1;
  x2 = y2;
  if (x1 > x2) {
      Ref<X> tmp = x1;
      x1 = x2;
      x2 = tmp;
      y1 << x1;
      y2 << x2;
    }

  cout << "x1 == x1: " << (x1 == x1) << "(1)";
  cout << " x1 == x2: " << (x1 == x2) << "(0)";
  cout << " x1 == y1: " << (x1 == y1) << "(1)";
  cout << " x1 == y2: " << (x1 == y2) << "(0)" << endl;
  cout << "x1 != x1: " << (x1 != x1) << "(0)";
  cout << " x1 != x2: " << (x1 != x2) << "(1)";
  cout << " x1 != y1: " << (x1 != y1) << "(0)";
  cout << " x1 != y2: " << (x1 != y2) << "(1)" << endl;
  cout << "x1 <  x1: " << (x1 <  x1) << "(0)";
  cout << " x1 <  x2: " << (x1 <  x2) << "(1)";
  cout << " x1 <  y1: " << (x1 <  y1) << "(0)";
  cout << " x1 <  y2: " << (x1 <  y2) << "(1)" << endl;
  cout << "x1 <= x1: " << (x1 <= x1) << "(1)";
  cout << " x1 <= x2: " << (x1 <= x2) << "(1)";
  cout << " x1 <= y1: " << (x1 <= y1) << "(1)";
  cout << " x1 <= y2: " << (x1 <= y2) << "(1)" << endl;
  cout << "x1 >  x1: " << (x1 >  x1) << "(0)";
  cout << " x1 >  x2: " << (x1 >  x2) << "(0)";
  cout << " x1 >  y1: " << (x1 >  y1) << "(0)";
  cout << " x1 >  y2: " << (x1 >  y2) << "(0)" << endl;
  cout << "x1 >= x1: " << (x1 >= x1) << "(1)";
  cout << " x1 >= x2: " << (x1 >= x2) << "(0)";
  cout << " x1 >= y1: " << (x1 >= y1) << "(1)";
  cout << " x1 >= y2: " << (x1 >= y2) << "(0)" << endl;
}

void
test3()
{
   cout << "nx = " << X::nx << "(0) (outer scope on entry)" << endl;
  {
  Ref<X> x1 = new X;
  Ref<X> x2 = new X;
  X *x2p = x2.pointer();
  cout << "nx = " << X::nx << "(2) (inner scope after alloc)" << endl;
#if REF_MANAGE
  x2->unmanage();
#endif

  x2 = x1.pointer();
  if (x1 != x1) abort();
  if (x2 != x1) abort();
#if REF_MANAGE
  // x2 was unmanaged, so delete it manually
  delete x2p;
#endif
  cout << "nx = " << X::nx << "(1) (inner scope after assign)" << endl;
  

  int i;
  for (i=1000000; i; i--) {
      x1->reference();
    }
  for (i=1000000; i; i--) {
      x1->dereference();
    }
  for (i=1000000; i; i--) {
      Ref<X> y = x1;
    }

  cout << "nx = " << X::nx << "(1) (inner scope)" << endl;
  }

   cout << "nx = " << X::nx << "(0) (outer scope on exit)" << endl;
}

#if HAVE_PTHREAD==1 && REF_USE_LOCKS==1
static Ref<X> sx1, sx2;
void *
test4_run(void *)
{
  for (int i=0; i<20000; i++) {
      Ref<X> x1 = sx1;
      Ref<X> x2 = sx2;
      x1 = x2;
      x1 = x2;
      x1 = x1;
      x1 = 0;
    }
  return 0;
}
void
test4()
{
  cout << "test4: running pthread tests" << endl;
  sx1 = new X;
  sx2 = new X;
  sx1->use_locks(true);
  sx2->use_locks(true);
  for (int iter=0; iter<10; iter++) {
      const int nthread = 2;
      pthread_t id[4];
      for (int i=0; i<nthread; i++) {
          pthread_create(&id[i],0,test4_run,0);
        }
      for (int i=0; i<nthread; i++) {
          void *v;
          pthread_join(id[i],&v);
        }
    }
  sx1 = 0;
  sx2 = 0;
}
#else
void test4() {}
#endif

int
main()
{
  cout << "X::nx before test1: " << X::nx << "(0)" << endl;
  test1();
  cout << "X::nx before test2: " << X::nx << "(0)" << endl;
  test2();
  cout << "X::nx before test3: " << X::nx << "(0)" << endl;
  test3();
  cout << "X::nx before test4: " << X::nx << "(0)" << endl;
  test4();

  cout << "X::nx before exit: " << X::nx << "(0)" << endl;

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
