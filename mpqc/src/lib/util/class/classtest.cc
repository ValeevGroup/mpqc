//
// classtest.cc
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

// a simple program to test the class stuff

#include <iostream.h>

#include <util/class/class.h>

#undef SIMPLE_TEST

#define A_parents virtual_base public DescribedClass
class A: A_parents {
#define CLASSNAME A
#include <util/class/classd.h>
  private:
    int i;
  public:
    A():i(1) {};
    ~A() { cout << "A dtor\n"; };
};

#define CLASSNAME A
#define PARENTS A_parents
#include <util/class/classi.h>
void *
A::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  DescribedClass::_castdown(cd) ;
  return do_castdowns(casts,cd);
}

#ifndef SIMPLE_TEST

#define B_parents public A
class B: B_parents {
#define CLASSNAME B 
#include <util/class/classd.h>
  private:
    int ib;
  public:
    B():ib(2) {};
    ~B() { cout << "B dtor\n"; };
};

#define CLASSNAME B
#define PARENTS B_parents
#include <util/class/classi.h>
void *
B::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =   A::_castdown(cd) ;
  return do_castdowns(casts,cd);
}

#define C_parents virtual_base public DescribedClass
class C: C_parents {
#define CLASSNAME C 
#include <util/class/classd.h>
  private:
    int i;
  public:
    C():i(3) {};
    ~C() { cout << "C dtor\n"; };
};

#define CLASSNAME C
#define PARENTS C_parents
#include <util/class/classi.h>
void *
C::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =   DescribedClass::_castdown(cd) ;
  return do_castdowns(casts,cd);
}

#ifndef NO_VIRTUAL_BASES
#define D_parents public B, public C
class D: D_parents {
#define CLASSNAME D
#define HAVE_CTOR
#include <util/class/classd.h>
  private:
    int id;
    A* atst;
  public:
    D():id(4),atst(new A) {};
    ~D() { delete atst; cout << "D dtor\n"; };
};

#define CLASSNAME D
#define PARENTS D_parents
#define HAVE_CTOR
#include <util/class/classi.h>
void *
D::_castdown(const ClassDesc*cd)
{
  void* casts[2];
  casts[0] =  B::_castdown(cd);
  casts[1] =  C::_castdown(cd);
  return do_castdowns(casts,cd);
}
#endif

#endif /* ! SIMPLE_TEST */

DescribedClass_REF_dec(A);
#ifndef SIMPLE_TEST
DescribedClass_REF_dec(B);
DescribedClass_REF_dec(C);
#ifndef NO_VIRTUAL_BASES
DescribedClass_REF_dec(D);
#endif
#endif /* ! SIMPLE_TEST */

DescribedClass_REF_def(A);
#ifndef SIMPLE_TEST
DescribedClass_REF_def(B);
DescribedClass_REF_def(C);
#ifndef NO_VIRTUAL_BASES
DescribedClass_REF_def(D);
#endif
#endif /* ! SIMPLE_TEST */

main()
{
  ClassDesc::list_all_classes();

  cout << "using 0" << endl;
  const RefDescribedClass descl2(0);
  RefA aaa;
  cout << "getting aaaa" << endl;
  A* aaaa = 0; //aaa.pointer();
  cout << "using aaaa" << endl;
  const RefDescribedClass descl((aaaa==(A*)0)?(DescribedClass*)0:aaaa);
  cout << "using aaa.pointer()" << endl;
  const RefDescribedClass descl3((aaa.pointer()==(A*)0)?(DescribedClass*)0:aaa.pointer());

  A a;
  cout << "A name:" << a.class_name() << '\n';

#ifndef NO_VIRTUAL_BASES
  D* dtst = D::castdown(ClassDesc::name_to_class_desc("D")->create());

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

#ifndef SIMPLE_TEST
  D d;
  cout << "D name:" << d.class_name() << '\n';

  cout << "&d = " << (void*) &d << '\n';
  cout << "D::castdown(&d) = " << (void*) D::castdown(&d) << '\n';
  cout << "B::castdown(&d) = " << (void*) B::castdown(&d) << '\n';
  cout << "A::castdown(&d) = " << (void*) A::castdown(&d) << '\n';
  cout << "C::castdown(&d) = " << (void*) C::castdown(&d) << '\n';
  cout << "DescribedClass::castdown(&d) = "
       << (void*) DescribedClass::castdown(&d) << '\n';

  RefD dref(new D);
  RefA aref(dref);

  cout << "aref.pointer() is " << aref.pointer() << '\n';
  cout << "dref.pointer() is " << dref.pointer() << '\n';
  cout << "aref == dref gives " << (aref == dref) << '\n';

  dref.operator=(aref);

  cout << "aref.pointer() is " << aref.pointer() << '\n';
  cout << "dref.pointer() is " << dref.pointer() << '\n';
  cout << "aref == dref gives " << (aref == dref) << '\n';
#endif /* ! SIMPLE_TEST */

#endif // ! NO_VIRTUAL_BASES

}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
