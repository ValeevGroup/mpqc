//
// keyvaltest.cc
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

#include <iostream.h>
#include <fstream.h>

#include <util/class/class.h>
#include <util/keyval/ipv2.h>
#include <util/keyval/keyval.h>

#define A_parents virtual public DescribedClass
class A: A_parents {
#define CLASSNAME A
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classd.h>
  private:
    int i;
  public:
    A();
    A(const RefKeyVal&keyval);
    inline int& a() { return i; };
    virtual void print (ostream&s = cout)
    {
      s << "A::a = " << a() << '\n';
    }
};
DescribedClass_REF_dec(A);
DescribedClass_REF_def(A);
A::A():
  i(1)
{
}
A::A(const RefKeyVal& keyval):
  i(keyval->intvalue("a"))
{
}

#define CLASSNAME A
#define PARENTS A_parents
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
A::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =   DescribedClass::_castdown(cd) ;
  return do_castdowns(casts,cd);
}

#define B_parents public A
class B: B_parents {
#define CLASSNAME B 
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classd.h>
  private:
    int b_;
  public:
    B();
    B(const RefKeyVal&keyval);
    inline int& b() { return b_; };
    virtual void print (ostream&s = cout)
    {
      A::print(s);
      s << "B::b = " << b() << '\n';
    }
};
DescribedClass_REF_dec(B);
DescribedClass_REF_def(B);
B::B():
  b_(2)
{
}
B::B(const RefKeyVal&keyval):
  A(new PrefixKeyVal(keyval,"A")),
  b_(keyval->intvalue("b"))
{
}

#define CLASSNAME B
#define PARENTS B_parents
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
B::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  A::_castdown(cd) ;
  return do_castdowns(casts,cd);
}

#define C_parents virtual public DescribedClass
class C: C_parents {
#define CLASSNAME C 
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classd.h>
  private:
    int i;
  public:
    C();
    C(const RefKeyVal&keyval);
    inline int& c() { return i; };
    virtual void print (ostream&s = cout)
    {
      s << "C::c = " << c() << '\n';
    }
};
DescribedClass_REF_dec(C);
DescribedClass_REF_def(C);
C::C():
  i(3)
{
}
C::C(const RefKeyVal&keyval):
  i(keyval->intvalue("c"))
{
}

#define CLASSNAME C
#define PARENTS C_parents
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
C::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =   DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

#define D_parents public B, public C
class D: D_parents {
#define CLASSNAME D
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classd.h>
  private:
    int d_;
    RefA d_a_;
    RefB d_b_;
  public:
    D();
    D(const RefKeyVal&keyval);
    inline int& d() { return d_; }
    inline RefA da() { return d_a_; }
    inline RefB db() { return d_b_; }
    virtual void print (ostream&s = cout)
    {
      B::print(s);
      C::print(s);
      s << "D::a:\n";  da()->print(s);
      if ( (A*)d_a_.pointer() == A::castdown(db().pointer()) ) {
          cout << "a == b\n";
        }
      else {
          s << "D::b:\n";  db()->print(s);
        }
      s << "D::d = " << d() << '\n';
    }
};
DescribedClass_REF_dec(D);
DescribedClass_REF_def(D);
D::D():
  d_(4)
{
}
D::D(const RefKeyVal&keyval):
  B(new PrefixKeyVal(keyval,"B")),
  C(new PrefixKeyVal(keyval,"C")),
  d_(keyval->intvalue("d")),
  d_a_(A::castdown(keyval->describedclassvalue("a"))),
  d_b_(B::castdown(keyval->describedclassvalue("b")))
{
}

#define CLASSNAME D
#define PARENTS D_parents
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
D::_castdown(const ClassDesc*cd)
{
  void* casts[2];
  casts[0] =  B::_castdown(cd);
  casts[1] =  C::_castdown(cd);
  return do_castdowns(casts,cd);
}

main()
{
  ClassDesc::list_all_classes();

  // test IPV2
  IPV2::Status err;
  ifstream in(SRCDIR "/keyvaltest.in",ios::in);
  if (in.bad()) {
      cout << "couldn't open " << SRCDIR << "/keyvaltest.in" << endl;
      abort();
    }
  IPV2 *ipv2 = new IPV2();
  ipv2->read(in,cout);
  ipv2->print_tree(cout);
  const char* test = 0;
  ipv2->value_v(":forref:nest:x",&test,0,0);
  cout << "test = \"" << test << "\"" << endl;
  err = ipv2->truekeyword_v(":forref:a",&test,0,0);
  cout << "test = \"" << test << "\" (" << ipv2->error_message(err) << ")"
       << endl;
  err = ipv2->truekeyword_v(":forref:nest:x",&test,0,0);
  cout << "test = \"" << test << "\" (" << ipv2->error_message(err) << ")"
       << endl;
  err = ipv2->truekeyword_v(":forref:x",&test,0,0);
  cout << "test = \"" << test << "\" (" << ipv2->error_message(err) << ")"
       << endl;
  delete ipv2;
  ipv2 = 0;

  // test the test classes

  A a;
  cout << "A name:" << a.class_name() << '\n';

  D d;
  cout << "D name:" << d.class_name() << '\n';

  cout << "&d = " << (void*) &d << '\n';
  cout << "D::castdown(&d) = " << (void*) D::castdown(&d) << '\n';
  cout << "B::castdown(&d) = " << (void*) B::castdown(&d) << '\n';
  cout << "A::castdown(&d) = " << (void*) A::castdown(&d) << '\n';
  cout << "C::castdown(&d) = " << (void*) C::castdown(&d) << '\n';
  cout << "DescribedClass::castdown(&d) = "
       << (void*) DescribedClass::castdown(&d) << '\n';


  RefAssignedKeyVal akv = new AssignedKeyVal;

  akv->assign(":x",1);
  akv->assign(":y",3.0);

#define stringize(arg) # arg
#define show( arg ) do{cout<<"   " stringize(arg) "="<<(arg);}while(0)

  show( akv->exists(":x") );  show( akv->errormsg() ); cout << '\n';
  show( akv->exists(":z") );  show (akv->errormsg() ); cout << '\n';
  show( akv->intvalue(":y") );  show( akv->errormsg() ); cout << '\n';
  show( akv->doublevalue(":x") );  show( akv->errormsg() ); cout << '\n';
  show( akv->intvalue(":x") );  show (akv->errormsg() ); cout << '\n';
  show( akv->intvalue("x") );  show (akv->errormsg() ); cout << '\n';
  show( akv->intvalue(":z") );  show (akv->errormsg() ); cout << '\n';

  RefKeyVal pkv = new ParsedKeyVal(SRCDIR "/keyvaltest.in");

  show( pkv->exists(":x") );  show( pkv->errormsg() ); cout << '\n';
  show( pkv->exists(":z") );  show (pkv->errormsg() ); cout << '\n';
  show( pkv->intvalue(":y") );  show( pkv->errormsg() ); cout << '\n';
  show( pkv->doublevalue(":x") );  show( pkv->errormsg() ); cout << '\n';
  show( pkv->intvalue(":x") );  show (pkv->errormsg() ); cout << '\n';
  show( pkv->intvalue("x") );  show (pkv->errormsg() ); cout << '\n';
  show( pkv->intvalue(":z") );  show (pkv->errormsg() ); cout << '\n';

  show ( pkv->exists("test:object_d") ); show(pkv->errormsg()); cout << '\n';

  RefDescribedClass rdc = pkv->describedclassvalue("test:object");
  show (pkv->errormsg() ); cout << '\n';
  show( rdc.pointer() ); cout << '\n';
  RefA ra(rdc);
  show( ra.pointer() ); cout << '\n';

  show( pkv->intvalue(":test:object:d") ); cout << '\n';

  pkv->dump();

  show( ra.pointer() ); cout << '\n';
  if (ra.nonnull()) { ra->print(); cout << '\n'; }

  cout << "Testing string keyvals" << endl;
  RefParsedKeyVal strkv = new ParsedKeyVal();
  strkv->parse_string("<B>:(b=123456)");
  RefB strb(strkv->describedclassvalue());
  if (strb.nonnull()) { strb->print(); cout << endl; }

  strkv = new ParsedKeyVal();
  strkv->parse_string("TOP=24");
  int strkvint = strkv->intvalue();
  cout << "strkvint = " << strkvint << endl;

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
