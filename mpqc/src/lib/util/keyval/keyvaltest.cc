
// a simple program to test the class stuff

#include <iostream.h>

#include <util/class/class.h>
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
    A(KeyVal&keyval);
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
A::A(KeyVal&keyval):
  i(keyval.intvalue("a"))
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
  void* casts[] =  { DescribedClass::_castdown(cd) };
  return do_castdowns(casts,cd);
}

#define B_parents public A
class B: B_parents {
#define CLASSNAME B 
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classd.h>
  private:
    int i;
  public:
    B();
    B(KeyVal&keyval);
    inline int& b() { return i; };
    virtual void print (ostream&s = cout)
    {
      A::print(s);
      s << "B::b = " << b() << '\n';
    }
};
DescribedClass_REF_dec(B);
DescribedClass_REF_def(B);
B::B():
  i(2)
{
}
B::B(KeyVal&keyval):
  A(PrefixKeyVal("A",keyval)),
  i(keyval.intvalue("b"))
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
  void* casts[] =  { A::_castdown(cd) };
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
    C(KeyVal&keyval);
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
C::C(KeyVal&keyval):
  i(keyval.intvalue("c"))
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
  void* casts[] =  { DescribedClass::_castdown(cd) };
  return do_castdowns(casts,cd);
}

#define D_parents public B, public C
class D: D_parents {
#define CLASSNAME D
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classd.h>
  private:
    int i;
    RefA _a;
    RefB _b;
  public:
    D();
    D(KeyVal&keyval);
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
DescribedClass_REF_dec(D);
DescribedClass_REF_def(D);
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

#define CLASSNAME D
#define PARENTS D_parents
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
D::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { B::_castdown(cd), C::_castdown(cd) };
  return do_castdowns(casts,cd);
}

main()
{
  ClassDesc::list_all_classes();

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

  ParsedKeyVal pkv("keyvaltest.in");

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
  RefA ra(rdc);
  show( ra.pointer() ); cout << '\n';

  show( pkv.intvalue(":test:object:d") ); cout << '\n';

  pkv.dump();

  show( ra.pointer() ); cout << '\n';
  if (ra.nonnull()) { ra->print(); cout << '\n'; }

}
