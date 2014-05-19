//
// classtest.cc
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

// a simple program to test the class stuff

#include <iostream>

#include <util/misc/formio.h>
#include <util/class/class.h>
#include <util/misc/scexception.h>

using namespace std;
using namespace sc;

#undef SIMPLE_TEST

class A: virtual public DescribedClass {
  private:
    int i;
  public:
    A():i(1) {};
    ~A() { cout << "A dtor\n"; };
};

static ClassDesc A_cd(typeid(A),"A",1,"virtual public DescribedClass");

#ifndef SIMPLE_TEST

class B: public A {
  private:
    int ib;
  public:
    B():ib(2) {};
    ~B() { cout << "B dtor\n"; };
};

static ClassDesc B_cd(typeid(B),"B",1,"public A");

class C: virtual public DescribedClass {
  private:
    int i;
  public:
    C():i(3) {};
    ~C() { cout << "C dtor\n"; };
};

static ClassDesc C_cd(typeid(C),"C",1,"virtual public DescribedClass");

class D: public B, public C {
  private:
    int id;
    A* atst;
  public:
    D():id(4),atst(new A) {};
    ~D() { delete atst; cout << "D dtor\n"; };
};

static ClassDesc D_cd(typeid(D),"D",1,"public B, public C",create<D>);

#endif /* ! SIMPLE_TEST */

int main(int argc, char* argv[])
{
  ClassDesc::list_all_classes();

  // try to construct duplicate ClassDesc for A
  try {
    ClassDesc A_cd(typeid(A),"A",1,"virtual public DescribedClass");
  }
  catch(sc::ProgrammingError& e) {
    cout << "tried constructing duplicate ClassDesc for A, caught ProgrammingError (as expected):" << endl;
    cout << e.what() << endl;
  }

  // try to construct duplicate ClassDesc for struct A but in different scope using same name
  try {
    struct A {};
    ClassDesc cd(typeid(A),"A",1,"");
  }
  catch(sc::ProgrammingError& e) {
    cout << "tried constructing ClassDesc for another class A in a different scope, caught ProgrammingError (as expected):" << endl;
    cout << e.what() << endl;
  }

  cout << indent << "using 0" << endl;
  const Ref<DescribedClass> descl2(0);
  Ref<A> aaa;
  cout << "getting aaaa" << endl;
  A* aaaa = 0; //aaa.pointer();
  cout << "using aaaa" << endl;
  const Ref<DescribedClass> descl((aaaa==(A*)0)?(DescribedClass*)0:aaaa);
  cout << "using aaa.pointer()" << endl;
  const Ref<DescribedClass> descl3((aaa.pointer()==(A*)0)?(DescribedClass*)0:aaa.pointer());

  A a;
  cout << "A name:" << a.class_name() << '\n';

  D* dtst = dynamic_cast<D*>(ClassDesc::name_to_class_desc("D")->create());
  delete dtst;

  // check the compiler's handling of virtual inheritance
  D* dt = new D;
  C* ct = dt;
  B* bt = dt;
  cout << "virtual inheritance test:" << endl;
  dt->reference();
  cout << "The following three numbers should be equal:" << endl;
  cout << ' ' << dt->nreference()
       << ' ' << ct->nreference()
       << ' ' << bt->nreference() << endl;
  ct->reference();
  cout << "The following three numbers should be equal:" << endl;
  cout << ' ' << dt->nreference()
       << ' ' << ct->nreference()
       << ' ' << bt->nreference() << endl;
  bt->reference();
  cout << "The following three numbers should be equal:" << endl;
  cout << ' ' << dt->nreference()
       << ' ' << ct->nreference()
       << ' ' << bt->nreference() << endl;
  cout << "done with virtual inheritance test:" << endl;
  dt->dereference();
  if (dt->nreference() == 0) delete dt;
  ct->dereference();
  if (ct->nreference() == 0) delete ct;
  bt->dereference();
  if (bt->nreference() == 0) delete bt;

#ifndef SIMPLE_TEST
  D d;
  cout << "D name:" << d.class_name() << '\n';

  cout << "&d = " << (void*) &d << '\n';
  cout << "dynamic_cast<D*>(&d) = " << (void*) dynamic_cast<D*>(&d) << '\n';
  cout << "dynamic_cast<B*>(&d) = " << (void*) dynamic_cast<B*>(&d) << '\n';
  cout << "dynamic_cast<A*>(&d) = " << (void*) dynamic_cast<A*>(&d) << '\n';
  cout << "dynamic_cast<C*>(&d) = " << (void*) dynamic_cast<C*>(&d) << '\n';
  cout << "dynamic_cast<DescribedClass*>(&d) = "
       << (void*) dynamic_cast<DescribedClass*>(&d) << '\n';

  Ref<D> dref(new D);
  Ref<A> aref(dref);

  cout << "aref.pointer() is " << aref.pointer() << '\n';
  cout << "dref.pointer() is " << dref.pointer() << '\n';
  cout << "aref == dref gives " << (aref == dref) << '\n';

  dref << aref;

  cout << "aref.pointer() is " << aref.pointer() << '\n';
  cout << "dref.pointer() is " << dref.pointer() << '\n';
  cout << "aref == dref gives " << (aref == dref) << '\n';
#endif /* ! SIMPLE_TEST */

}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
