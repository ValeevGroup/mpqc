
// a simple program to test the state stuff

#include <iostream.h>

#include <util/class/class.h>
#include <util/keyval/keyval.h>
#include <util/state/state.h>

#define A_parents virtual public SavableState
class A: A_parents {
#   define CLASSNAME A
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    int i;
    int* array;
    double d;
  public:
    A();
    A(KeyVal&);
    A(StateIn&);
    void save_data_state(StateOut&);
    inline int& a() { return i; };
    virtual void print (ostream&s = cout)
    {
      s << "A::a = " << a() << '\n';
      s << "A::array = {"
        << array[0] << ' '
        << array[1] << ' '
        << array[2] << ' '
        << array[3]
        << "}\n";
    }
};
SavableState_REF_dec(A);
SavableState_REF_def(A);
A::A():
  i(1),
  array(new int[4]),
  d(-1.24)
{
  array[0] = 4;
  array[1] = 3;
  array[2] = 2;
  array[3] = 1;
}
A::A(KeyVal&keyval):
  i(keyval.intvalue("a")),
  array(new int[4]),
  d(-1.24)

{
  array[0] = 4;
  array[1] = 3;
  array[2] = 2;
  array[3] = 8;
}
A::A(StateIn&s):
  SavableState(s,class_desc_)
{
  char* junk;
  s.get(d);
  s.getstring(junk); delete[] junk;
  s.get(i);
  s.getstring(junk); delete[] junk;
  s.get(array);
}
void
A::save_data_state(StateOut&s)
{
  const char* t1 = "test string";
  const char* t2 = "test2\nstring";
  char* t1c = strcpy(new char[strlen(t1)+1],t1);
  char* t2c = strcpy(new char[strlen(t2)+1],t2);
  s.put(d);
  s.putstring(t1c);
  s.put(i);
  s.putstring(t2c);
  s.put(array,4);
  delete[] t1c;
  delete[] t2c;
}

#define CLASSNAME A
#define PARENTS A_parents
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
A::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}

#define B_parents public A
class B: B_parents {
#   define CLASSNAME B 
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    int i;
  public:
    B();
    B(KeyVal&);
    B(StateIn&);
    void save_data_state(StateOut&);
    inline int& b() { return i; };
    virtual void print (ostream&s = cout)
    {
      A::print(s);
      s << "B::b = " << b() << '\n';
    }
};
SavableState_REF_dec(B);
SavableState_REF_def(B);
B::B():
  i(2)
{
}
B::B(KeyVal&keyval):
  A(PrefixKeyVal("A",keyval)),
  i(keyval.intvalue("b"))
{
}
B::B(StateIn&s):
  SavableState(s,class_desc_),
  A(s)
{
  s.get(i);
}
void
B::save_data_state(StateOut&s)
{
  A::save_data_state(s);
  s.put(i);
}

#define CLASSNAME B
#define PARENTS B_parents
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
B::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { A::_castdown(cd) };
  return do_castdowns(casts,cd);
}

#define C_parents virtual public SavableState
class C: C_parents {
#   define CLASSNAME C 
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    int i;
  public:
    C();
    C(KeyVal&keyval);
    C(StateIn&);
    void save_data_state(StateOut&);
    inline int& c() { return i; };
    virtual void print (ostream&s = cout)
    {
      s << "C::c = " << c() << '\n';
    }
};
SavableState_REF_dec(C);
SavableState_REF_def(C);
C::C():
  i(3)
{
}
C::C(KeyVal&keyval):
  i(keyval.intvalue("c"))
{
}
C::C(StateIn&s):
  SavableState(s,class_desc_)
{
  s.get(i);
}
void
C::save_data_state(StateOut&s)
{
  s.put(i);
}

#define CLASSNAME C
#define PARENTS C_parents
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
C::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}

#define D_parents public B, public C
class D: D_parents {
#   define CLASSNAME D
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    int i;
    RefA _a;
    RefB _b;
  public:
    D();
    D(KeyVal&);
    D(StateIn&);
    void save_data_state(StateOut&);
    inline int& d() { return i; }
    inline RefA da() { return _a; }
    inline RefB db() { return _b; }
    virtual void print (ostream&s = cout)
    {
      B::print(s);
      C::print(s);
      s << "D::a:\n";  da()->print(s);
      if ( _a == A::castdown(db().pointer()) ) {
          cout << "a == b\n";
        }
      else {
          s << "D::b:\n";  db()->print(s);
        }
      s << "D::d = " << d() << '\n';
    }
};
SavableState_REF_dec(D);
SavableState_REF_def(D);
D::D():
  i(4)
{
}
D::D(KeyVal&keyval):
  B(PrefixKeyVal("B",keyval)),
  C(PrefixKeyVal("C",keyval)),
  i(keyval.intvalue("d")),
  _a(A::castdown(keyval.describedclassvalue("a"))),
  _b(B::castdown(keyval.describedclassvalue("b")))
{
}
D::D(StateIn&s):
  SavableState(s,class_desc_),
  B(s),
  C(s)
{
  int* nullref;
  s.get(nullref);
  s.get(i);
  _a = A::restore_state(s);
  _b = B::restore_state(s);
}
void
D::save_data_state(StateOut&s)
{
  B::save_data_state(s);
  C::save_data_state(s);
  s.put((int*)0,0);
  s.put(i);
  _a->save_state(s);
  _b->save_state(s);
}

#define CLASSNAME D
#define PARENTS D_parents
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
D::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { B::_castdown(cd), C::_castdown(cd) };
  return do_castdowns(casts,cd);
}

#define StateOutType StateOutText
#define StateInType StateInText

//#define StateOutType StateOutBinXDR
//#define StateInType StateInBinXDR

main()
{
  RefA ra;

  ClassDesc::list_all_classes();

  ra = 0;

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

  AssignedKeyVal akv;

  akv.assign(":x",1);
  akv.assign(":y",3.0);

#define stringize(arg) # arg
#define show( arg ) do{cout<<"   " stringize(arg) "="<<(arg);}while(0)

  show( akv.exists(":x") );  show( akv.errormsg() ); cout << '\n';
  show( akv.exists(":z") );  show (akv.errormsg() ); cout << '\n';
  show( akv.intvalue(":y") );  show( akv.errormsg() ); cout << '\n';
  show( akv.doublevalue(":x") );  show( akv.errormsg() ); cout << '\n';
  show( akv.intvalue(":x") );  show (akv.errormsg() ); cout << '\n';
  show( akv.intvalue("x") );  show (akv.errormsg() ); cout << '\n';
  show( akv.intvalue(":z") );  show (akv.errormsg() ); cout << '\n';

  ParsedKeyVal pkv("statetest.in");

  show( pkv.exists(":x") );  show( pkv.errormsg() ); cout << '\n';
  show( pkv.exists(":z") );  show (pkv.errormsg() ); cout << '\n';
  show( pkv.intvalue(":y") );  show( pkv.errormsg() ); cout << '\n';
  show( pkv.doublevalue(":x") );  show( pkv.errormsg() ); cout << '\n';
  show( pkv.intvalue(":x") );  show (pkv.errormsg() ); cout << '\n';
  show( pkv.intvalue("x") );  show (pkv.errormsg() ); cout << '\n';
  show( pkv.intvalue(":z") );  show (pkv.errormsg() ); cout << '\n';

  RefDescribedClass rdc = pkv.describedclassvalue("test:object");
  show (pkv.errormsg() ); cout << '\n';
  show( rdc.pointer() ); cout << '\n';
  ra = A::castdown(rdc);
  show( ra.pointer() ); cout << '\n';

  show( pkv.intvalue(":test:object:d") ); cout << '\n';

  //pkv.dump();

  show( ra.pointer() ); cout << '\n';
  if (ra.nonnull()) { ra->print(); cout << '\n'; }

  ////////////////////////////////////////////////////////////////////
  // state tests

  cout << " -- saving state --\n";

  StateOutType soa("statetest.a.out");
  ra = new A(PrefixKeyVal("test:object_a",pkv));
  ra->save_object_state(soa);
  soa.forget();
  ra->save_object_state(soa);
  soa.flush();
  ra = A::castdown(rdc);
  StateOutText so("statetest.out");
  ra.save_state(so);
  RefA ra2;
  ra2.save_state(so);
  so.flush();

  cout << " -- restoring state --\n";

  StateInType sia("statetest.a.out");
  ra = new A(sia);
  ra = new A(sia);
  if (ra.nonnull()) { ra->print(); cout << '\n'; }
  StateInText si("statetest.out");
  //ra = A::restore_state(si);
  ra.restore_state(si);
  ra2.restore_state(si);
  if (ra.nonnull()) { ra->print(); cout << '\n'; }
  cout << "ra2.nonnull() = " << ra2.nonnull() << "(should be 0)\n";
}
