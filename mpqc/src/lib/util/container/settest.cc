//
// settest.cc
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

#include <iomanip.h>

#include <util/ref/ref.h>
#include <util/container/array.h>
#include <util/container/set.h>

class Intg: public VRefCount {
 private:
  int _int;
 public:
  Intg();
  Intg(int i);
  void print();
};
Intg::Intg():_int(0) {};
Intg::Intg(int i):
_int(i)
{
};
void Intg::print() { cout << _int << endl; }

REF_dec(Intg);
REF_def(Intg);

ARRAY_dec(RefIntg);
ARRAY_def(RefIntg);

SET_dec(RefIntg);
SET_def(RefIntg);

ARRAYSET_dec(RefIntg);
ARRAYSET_def(RefIntg);

class A: public VRefCount {
    int a;
  public:
    virtual ~A() {};
};
REF_dec(A);
REF_def(A);
class B: public A {
    int b;
};
REF_dec(B);
REF_def(B);

main()
{
  RefA aar(new A);
  aar = 0;

  A* aap = new A;
  cout << "0x" << setbase(16) << setw(8) << setfill('0') << aap
       << "ref count = " << setbase(10) << aap->nreference() << endl;
//   delete aap;
//   cout << "0x" << setbase(16) << setw(8) << setfill('0') << aap
//        << "ref count = " << setbase(10) << aap->_reference_count_ << endl;
  aar = aap;
//   delete aap;
//   cout << "0x" << setbase(16) << setw(8) << setfill('0') << aap
//        << "ref count = " << setbase(10) << aap->_reference_count_ << endl;
  aar = 0;
  
  RefB b;
  //RefA a(b); // illegal
  RefA a(b.pointer());
  
  Intg* i1 = new Intg(101);
  i1->print();

  RefIntg ii = new Intg(100);

  ArraysetRefIntg as;

  int i;
  for (i=0; i<2; i++) {
      as.add(new Intg(i));
    }
  as[1] = ii;
  for (i=0; i<2; i++) {
      as[i]->print();
    }
  ii->print();
  i1->print();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
