//
// statetest.cc
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

// a simple program to test the state stuff

#include <iostream>

#include <util/misc/formio.h>

#include <util/class/class.h>
#include <util/keyval/keyval.h>
#include <util/state/state.h>

#include <util/state/linkage.h>

using namespace std;
using namespace sc;

#if 0 // normally 0
#  define StateOutTypeA StateOutText
#  define StateInTypeA StateInText
#  include <util/state/state_text.h>
#else
#  define StateOutTypeA StateOutBin
#  define StateInTypeA StateInBin
#  include <util/state/state_bin.h>
#endif

#if 0 // normally 0
#  define StateOutTypeB StateOutBin
#  define StateInTypeB StateInBin
#  include <util/state/state_bin.h>
#else
#  define StateOutTypeB StateOutText
#  define StateInTypeB StateInText
#  include <util/state/state_text.h>
#endif

class A: virtual public SavableState {
  private:
    int ia;
    int* array;
    double d;
    char *t1c;
    char *t2c;
  public:
    A();
    A(const Ref<KeyVal>&);
    A(StateIn&);
    ~A();
    void save_data_state(StateOut&);
    inline int& a() { return ia; };
    virtual void print (ostream&s = ExEnv::out0())
    {
      s << "A::t1c = " << t1c << '\n';
      s << "A::t2c = " << t2c << '\n';
      s << "A::a = " << a() << '\n';
      s << "A::d = " << d << '\n';
      s << "A::array = {"
        << array[0] << ' '
        << array[1] << ' '
        << array[2] << ' '
        << array[3]
        << "}\n";
    }
};

A::A():
  ia(1),
  array(new int[4]),
  d(-1.24)
{
  array[0] = 4;
  array[1] = 3;
  array[2] = 2;
  array[3] = 1;
  const char* t1 = "test string";
  const char* t2 = "test2\nstring";
  t1c = strcpy(new char[strlen(t1)+1],t1);
  t2c = strcpy(new char[strlen(t2)+1],t2);
}
A::A(const Ref<KeyVal>&keyval):
  ia(keyval->intvalue("a")),
  array(new int[4]),
  d(-1.24)

{
  array[0] = 4;
  array[1] = 3;
  array[2] = 2;
  array[3] = 8;
  const char* t1 = "test string";
  const char* t2 = "test2\nstring";
  t1c = strcpy(new char[strlen(t1)+1],t1);
  t2c = strcpy(new char[strlen(t2)+1],t2);
}
A::A(StateIn&s):
  SavableState(s)
{
  s.get(d,"d");
  s.getstring(t1c);
  s.get(ia,"a");
  s.getstring(t2c);
  s.get(array);
}
A::~A()
{
  delete[] array;
  delete[] t1c;
  delete[] t2c;
}
void
A::save_data_state(StateOut&s)
{
  s.put(d);
  s.putstring(t1c);
  s.put(ia);
  s.putstring(t2c);
  s.put(array,4);
}

static ClassDesc A_cd(typeid(A),"A",1,"virtual public SavableState",
                      create<A>, create<A>, create<A>);

class B: public A {
  private:
    int ib;
  public:
    B();
    B(const Ref<KeyVal>&);
    B(StateIn&);
    void save_data_state(StateOut&);
    inline int& b() { return ib; };
    virtual void print (ostream&s = cout)
    {
      A::print(s);
      s << "B::b = " << b() << '\n';
    }
};

B::B():
  ib(2)
{
}
B::B(const Ref<KeyVal>&keyval):
  A(keyval),
  ib(keyval->intvalue("b"))
{
}
B::B(StateIn&s):
  SavableState(s),
  A(s)
{
  s.get(ib);
}
void
B::save_data_state(StateOut&s)
{
  A::save_data_state(s);
  s.put(ib);
}

static ClassDesc B_cd(typeid(B),"B",1,"public A",
                      create<B>,create<B>,create<B>);

class C: virtual public SavableState {
  private:
    int ic;
  public:
    C();
    C(const Ref<KeyVal>&keyval);
    C(StateIn&);
    void save_data_state(StateOut&);
    inline int& c() { return ic; };
    virtual void print (ostream&s = cout)
    {
      s << "C::c = " << c() << '\n';
    }
};

C::C():
  ic(3)
{
}
C::C(const Ref<KeyVal>&keyval):
  ic(keyval->intvalue("c"))
{
}
C::C(StateIn&s):
  SavableState(s)
{
  s.get(ic);
}
void
C::save_data_state(StateOut&s)
{
  s.put(ic);
}

static ClassDesc C_cd(typeid(C),"C",1,"virtual public SavableState",
                      create<C>,create<C>,create<C>);

class D: public B, public C {
  private:
    int id;
    long ld;
    char cd;
    float fd;
    double dd;
    Ref<A> _a;
    Ref<B> _b;
    char *cdat;
    int *idat;
    long *ldat;
    float *fdat;
    double *ddat;
    std::string sdat;
  public:
    D();
    D(const Ref<KeyVal>&);
    D(StateIn&);
    ~D();
    void save_data_state(StateOut&);
    inline int& d() { return id; }
    inline Ref<A> da() { return _a; }
    inline Ref<B> db() { return _b; }
    virtual void print (ostream&s = cout)
    {
      B::print(s);
      C::print(s);
      s << "D::a:\n";
      if (da()) {
          da()->print(s);
        }
      else {
          s << "null\n";
        }
      if ( _a.pointer() == dynamic_cast<A*>(db().pointer())) 
        {
          cout << "a == b\n";
        }
      else {
          s << "D::b:\n";  db()->print(s);
        }
      s << "D::d = " << d() << '\n';
      s << "D::sdat = " << sdat << std::endl;
    }
};

D::D()
{
  id = 4;
  ld = 14;
  cd = 'd';
  fd = 4.1;
  dd = 8.2;
  ddat = new double[4];
  fdat = new float[4];
  idat = new int[4];
  ldat = new long[4];
  cdat = new char[4];
  cdat[0]=(cdat[1]=(cdat[2]=(cdat[3]='a')+1)+1)+1;
  idat[0]=(idat[1]=(idat[2]=(idat[3]=1)+1)+1)+1;
  ldat[0]=(ldat[1]=(ldat[2]=(ldat[3]=1)+1)+1)+1;
  fdat[0]=(fdat[1]=(fdat[2]=(fdat[3]=1.0)+1)+1)+1;
  ddat[0]=(ddat[1]=(ddat[2]=(ddat[3]=1.0)+1)+1)+1;
  sdat = "Test of std::string";
}
D::D(const Ref<KeyVal>&keyval):
  B(keyval),
  C(keyval),
  id(keyval->intvalue("di")),
  ld(keyval->longvalue("dl")),
  cd(keyval->charvalue("dc")),
  fd(keyval->floatvalue("df")),
  dd(keyval->doublevalue("dd")),
  _a(dynamic_cast<A*>(keyval->describedclassvalue("da").pointer())),
  _b(dynamic_cast<B*>(keyval->describedclassvalue("db").pointer()))
{
  ddat = new double[4];
  fdat = new float[4];
  idat = new int[4];
  ldat = new long[4];
  cdat = new char[4];
  cdat[0]=(cdat[1]=(cdat[2]=(cdat[3]='a')+1)+1)+1;
  idat[0]=(idat[1]=(idat[2]=(idat[3]=1)+1)+1)+1;
  ldat[0]=(ldat[1]=(ldat[2]=(ldat[3]=1)+1)+1)+1;
  fdat[0]=(fdat[1]=(fdat[2]=(fdat[3]=1.0)+1)+1)+1;
  ddat[0]=(ddat[1]=(ddat[2]=(ddat[3]=1.0)+1)+1)+1;
  sdat = "Test of std::string";
}
D::D(StateIn&s):
  SavableState(s),
  B(s),
  C(s)
{
  s.get(id,"di");
  s.get(ld,"dl");
  s.get(cd,"dc");
  s.get(fd,"df");
  s.get(dd,"dd");
  char *junk;
  s.getstring(junk);
  delete[] junk;
  _a << SavableState::key_restore_state(s,"da");
  s.getstring(junk);
  delete[] junk;
  _b << SavableState::key_restore_state(s,"db");
  s.get(ddat);
  s.get(fdat);
  s.get(idat);
  s.get(ldat);
  s.get(cdat);
  s.get(sdat);
}
void
D::save_data_state(StateOut&s)
{
  B::save_data_state(s);
  C::save_data_state(s);
  s.put(id);
  s.put(ld);
  s.put(cd);
  s.put(fd);
  s.put(dd);
  s.putstring("here begins _a");
  SavableState::save_state(_a.pointer(), s);
  s.putstring("here begins _b");
  SavableState::save_state(_b.pointer(),s);
  s.put(ddat,4);
  s.put(fdat,4);
  s.put(idat,4);
  s.put(ldat,4);
  s.put(cdat,4);
  s.put(sdat);
}
D::~D()
{
  delete[] ddat;
  delete[] fdat;
  delete[] idat;
  delete[] ldat;
  delete[] cdat;
}

static ClassDesc D_cd(typeid(D),"D",1,"public B, public C",
                      create<D>, create<D>, create<D>);

int
main(int argc, char* argv[])
{
  Ref<A> ra;

  ClassDesc::list_all_classes();

  ra = 0;

  A a;
  cout << "A name:" << a.class_name() << endl;

  D d;
  cout << "D name:" << d.class_name() << endl;

  cout << "&d = " << (void*) &d << endl;
  cout << "dynamic_cast<D*>(&d) = " << (void*) dynamic_cast<D*>(&d) << endl;
  cout << "dynamic_cast<B*>(&d) = " << (void*) dynamic_cast<B*>(&d) << endl;
  cout << "dynamic_cast<A*>(&d) = " << (void*) dynamic_cast<A*>(&d) << endl;
  cout << "dynamic_cast<C*>(&d) = " << (void*) dynamic_cast<C*>(&d) << endl;
  cout << "dynamic_cast<DescribedClass*>(&d) = "
       << (void*) dynamic_cast<DescribedClass*>(&d) << endl;

  Ref<AssignedKeyVal> akv (new AssignedKeyVal);

  akv->assign(":x",1);
  akv->assign(":y",3.0);

#define stringize(arg) # arg
#define show( arg ) do{cout<<"   " stringize(arg) "="<<(arg);}while(0)

  show( akv->exists(":x") );  show( akv->errormsg() ); cout << endl;
  show( akv->exists(":z") );  show (akv->errormsg() ); cout << endl;
  show( akv->intvalue(":y") );  show( akv->errormsg() ); cout << endl;
  show( akv->doublevalue(":x") );  show( akv->errormsg() ); cout << endl;
  show( akv->intvalue(":x") );  show (akv->errormsg() ); cout << endl;
  show( akv->intvalue("x") );  show (akv->errormsg() ); cout << endl;
  show( akv->intvalue(":z") );  show (akv->errormsg() ); cout << endl;

  Ref<KeyVal> pkv = new ParsedKeyVal(SRCDIR "/statetest.in");

  show( pkv->exists(":x") );  show( pkv->errormsg() ); cout << endl;
  show( pkv->exists(":z") );  show (pkv->errormsg() ); cout << endl;
  show( pkv->intvalue(":y") );  show( pkv->errormsg() ); cout << endl;
  show( pkv->doublevalue(":x") );  show( pkv->errormsg() ); cout << endl;
  show( pkv->intvalue(":x") );  show (pkv->errormsg() ); cout << endl;
  show( pkv->intvalue("x") );  show (pkv->errormsg() ); cout << endl;
  show( pkv->intvalue(":z") );  show (pkv->errormsg() ); cout << endl;

  Ref<DescribedClass> rdc = pkv->describedclassvalue("test:object");
  show (pkv->errormsg() ); cout << endl;
  show( rdc.pointer() ); cout << endl;
  ra = dynamic_cast<A*>(rdc.pointer());
  show( ra.pointer() ); cout << endl;

  show( pkv->intvalue(":test:object:d") ); cout << endl;

  //pkv->dump();

  show( ra.pointer() ); cout << endl;
  if (ra) { ra->print(); cout << endl; }

  ////////////////////////////////////////////////////////////////////
  // state tests

  cout << " ------------- saving state ----------------" << endl;

  cout << " --- saving to A ---" << endl;
  StateOutTypeA soa("statetest.a.out");
  ra = new A(new PrefixKeyVal(pkv,"test:object_a"));
  cout << "  first a" << endl;
  ra->save_object_state(soa);
  soa.forget_references();
  cout << "  second a" << endl;
  ra->save_object_state(soa);
  ra = dynamic_cast<A*>(rdc.pointer());
  ra->save_state(soa);
  soa.flush();
  soa.close();
  cout << " --- saving to B ---" << endl;
  StateOutTypeB so("statetest.out");
  SavableState::save_state(ra.pointer(),so);
  Ref<A> ra2;
  SavableState::save_state(ra2.pointer(),so);
  so.close();

  cout << " ------------- restoring state ----------------" << endl;

  cout << " --- restoring from A ---" << endl;
  StateInTypeA sia("statetest.a.out");
  cout << "  first a" << endl;
  ra = new A(sia);
  cout << "  second a" << endl;
  ra = new A(sia);
  cout << "  last object" << endl;
  ra << SavableState::restore_state(sia);
  if (ra) { ra->print(); cout << endl; }
  if (sia.use_directory()) {
      cout << " --- restoring from A's directory ---" << endl;
      ra << SavableState::dir_restore_state(sia,"B:1");
      cout << "B:1 classname = " << ra->class_name() << endl;
    }
  sia.close();
  cout << " --- restoring from B ---" << endl;
  StateInTypeB si("statetest.out");
  //ra = A::restore_state(si);
  ra << SavableState::restore_state(si);
  ra2 << SavableState::restore_state(si);
  if (ra) { ra->print(); cout << endl; }
  cout << "ra2 = " << ra2 << "(should be 0)\n";
  si.close();

  if (sia.use_directory()) {
      sia.open("statetest.a.out");
      ExEnv::out0() << indent
           << " --- restoring from A's directory (2) ---" << endl;
      ra << SavableState::dir_restore_state(sia,"B:1");
      ExEnv::out0() << indent
           << "B:1 classname = " << ra->class_name() << endl;
      Ref<A> ra3;
      ra3 << SavableState::dir_restore_state(sia,"B:1");
      ExEnv::out0() << indent
           <<"first B:1: " << (void*) ra.pointer()
           << " second B:1: " << (void*) ra3.pointer()
           << endl;
    }
  ExEnv::out0() << indent << "objects in sia" << endl;
  sia.list_objects();

  if (sia.use_directory()) {
      cout << " ----- proxy tests ----- " << endl;
      Ref<D> d1; d1 << pkv->describedclassvalue("test2:proxy1");
      Ref<D> d2; d2 << pkv->describedclassvalue("test2:proxy2");
      cout << "d1 = " << (void*)d1.pointer()
           << " d2 = " << (void*)d2.pointer() << endl;
      if (d1) d1->print();
    }

  return 0;
}
