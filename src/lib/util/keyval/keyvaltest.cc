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

#include <iostream>
#include <fstream>

#include <util/misc/formio.h>
#include <util/class/class.h>
#include <util/keyval/ipv2.h>
#include <util/keyval/keyval.h>

using namespace std;
using namespace sc;

class A: virtual public DescribedClass {
  private:
    int i;
  public:
    A();
    A(const Ref<KeyVal>&keyval);
    inline const int& a() const { return i; };
    virtual void print (ostream&s = cout) const
    {
      s << "A::a = " << a() << '\n';
    }
};

A::A():
  i(1)
{
}
A::A(const Ref<KeyVal>& keyval):
  i(keyval->intvalue("a"))
{
}

static ClassDesc A_cd(typeid(A),"A",1,"virtual public DescribedClass",
                      create<A>,create<A>,0);
 
class B: public A {
  private:
    int b_;
  public:
    B();
    B(const Ref<KeyVal>&keyval);
    inline const int& b() const { return b_; };
    virtual void print (ostream&s = cout) const
    {
      A::print(s);
      s << "B::b = " << b() << '\n';
    }
};

B::B():
  b_(2)
{
}
B::B(const Ref<KeyVal>&keyval):
  A(new PrefixKeyVal(keyval,"A")),
  b_(keyval->intvalue("b"))
{
}

class ClassDesc B_cd(typeid(B),"B",1,"public A",create<B>,create<B>);

class C: virtual public DescribedClass {
  private:
    int i;
  public:
    C();
    C(const Ref<KeyVal>&keyval);
    inline const int& c() const { return i; };
    virtual void print (ostream&s = cout) const
    {
      s << "C::c = " << c() << '\n';
    }
};

C::C():
  i(3)
{
}
C::C(const Ref<KeyVal>&keyval):
  i(keyval->intvalue("c"))
{
}

static ClassDesc C_cd(typeid(C),"C",1,"virtual public DescribedClass",
                      create<C>, create<C>);

class D: public B, public C {
  private:
    int d_;
    Ref<A> d_a_;
    Ref<B> d_b_;
  public:
    D();
    D(const Ref<KeyVal>&keyval);
    inline const int& d() const { return d_; }
    inline Ref<A> da() const { return d_a_; }
    inline Ref<B> db() const { return d_b_; }
    virtual void print (ostream&s = cout) const
    {
      B::print(s);
      C::print(s);
      s << "D::a:\n";  da()->print(s);
      if ( (A*)d_a_.pointer() == dynamic_cast<A*>(db().pointer()) ) {
          cout << "a == b\n";
        }
      else {
          s << "D::b:\n";  db()->print(s);
        }
      s << "D::d = " << d() << '\n';
    }
};

D::D():
  d_(4)
{
}
D::D(const Ref<KeyVal>&keyval):
  B(new PrefixKeyVal(keyval,"B")),
  C(new PrefixKeyVal(keyval,"C")),
  d_(keyval->intvalue("d")),
  d_a_(dynamic_cast<A*>(keyval->describedclassvalue("a").pointer())),
  d_b_(dynamic_cast<B*>(keyval->describedclassvalue("b").pointer()))
{
}

static ClassDesc D_cd(typeid(D),"D",1,"public B, public C",
                      create<D>,create<D>);

int
main(int argc, char* argv[])
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
  cout << "dynamic_cast<D*>(&d) = " << (void*) dynamic_cast<D*>(&d) << '\n';
  cout << "dynamic_cast<B*>(&d) = " << (void*) dynamic_cast<B*>(&d) << '\n';
  cout << "dynamic_cast<A*>(&d) = " << (void*) dynamic_cast<A*>(&d) << '\n';
  cout << "dynamic_cast<C*>(&d) = " << (void*) dynamic_cast<C*>(&d) << '\n';
  cout << "dynamic_cast<DescribedClass*>(&d) = "
       << (void*) dynamic_cast<DescribedClass*>(&d) << '\n';


  Ref<AssignedKeyVal> akv = new AssignedKeyVal;

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

  Ref<KeyVal> pkv = new ParsedKeyVal(SRCDIR "/keyvaltest.in");
  pkv->verbose(1);

  cout << "Initial unseen keywords:" << endl;
  cout << incindent;
  pkv->print_unseen(cout);
  cout << decindent;
  cout << "done" << endl;

  // size tests
  for (int i=0; i < pkv->count("memory"); i++) {
      cout << "memory:" << i << " = " << pkv->sizevalue("memory",i) << endl;
    }

  show( pkv->exists(":s1") );  show( pkv->errormsg() ); cout << '\n';
  show( pkv->exists(":s2") );  show( pkv->errormsg() ); cout << '\n';
  show( pkv->stringvalue(":s1") );  show( pkv->errormsg() ); cout << '\n';
  show( pkv->stringvalue(":s2") );  show( pkv->errormsg() ); cout << '\n';
  char *tmp;
  show( tmp=pkv->pcharvalue(":s1") );  show( pkv->errormsg() ); cout << '\n';
  delete[] tmp;
  show( tmp=pkv->pcharvalue(":s2") );  show( pkv->errormsg() ); cout << '\n';
  delete[] tmp;

  show( pkv->exists(":x") );  show( pkv->errormsg() ); cout << '\n';
  show( pkv->exists(":z") );  show (pkv->errormsg() ); cout << '\n';
  show( pkv->intvalue(":y") );  show( pkv->errormsg() ); cout << '\n';
  show( pkv->doublevalue(":x") );  show( pkv->errormsg() ); cout << '\n';
  show( pkv->intvalue(":x") );  show (pkv->errormsg() ); cout << '\n';
  show( pkv->intvalue("x") );  show (pkv->errormsg() ); cout << '\n';
  show( pkv->intvalue(":z") );  show (pkv->errormsg() ); cout << '\n';

  show ( pkv->exists("test:object_d") ); show(pkv->errormsg()); cout << '\n';

  Ref<DescribedClass> rdc = pkv->describedclassvalue("test:object");
  show (pkv->errormsg() ); cout << '\n';
  show( rdc.pointer() ); cout << '\n';
  Ref<A> ra; ra << rdc;
  show( ra.pointer() ); cout << '\n';

  show( pkv->intvalue(":test:object:d") ); cout << '\n';

  pkv->dump();

  cout << "Final unseen keywords:" << endl;
  pkv->exists("forref");
  pkv->exists("testintco");
  cout << incindent;
  pkv->print_unseen(cout);
  cout << decindent;
  cout << "done" << endl;

  show( ra.pointer() ); cout << '\n';
  if (ra.nonnull()) { ra->print(); cout << '\n'; }

  cout << "Testing string keyvals" << endl;
  Ref<ParsedKeyVal> strkv = new ParsedKeyVal();
  cout << "  parsing" << endl;
  strkv->parse_string("<B>:(b=123456)");
  cout << "  reading" << endl;
  Ref<B> strb; strb << strkv->describedclassvalue();
  if (strb.nonnull()) {
      cout << "  printing" << endl;
      strb->print(); cout << endl;
    }

  cout << "Testing parsed keyvals TOP keyword" << endl;
  strkv = new ParsedKeyVal();
  cout << "  parsing" << endl;
  strkv->parse_string("TOP=24");
  cout << "  reading" << endl;
  int strkvint = strkv->intvalue();
  cout << "strkvint = " << strkvint << endl;

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
