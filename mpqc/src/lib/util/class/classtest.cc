
// a simple program to test the class stuff

#include <iostream.h>

#include <util/class/class.h>

#undef SIMPLE_TEST

#define A_parents virtual public DescribedClass
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

#define C_parents virtual public DescribedClass
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

#define GNU_BUG_VIRTUAL virtual
#define D_parents GNU_BUG_VIRTUAL public B, GNU_BUG_VIRTUAL public C
class D: D_parents {
#define CLASSNAME D
#include <util/class/classd.h>
  private:
    int id;
  public:
    D():id(4) {};
    ~D() { cout << "D dtor\n"; };
};

#define CLASSNAME D
#define PARENTS D_parents
#include <util/class/classi.h>
void *
D::_castdown(const ClassDesc*cd)
{
  void* casts[2];
  casts[0] =  B::_castdown(cd);
  casts[1] =  C::_castdown(cd);
  return do_castdowns(casts,cd);
}

#endif /* ! SIMPLE_TEST */

DescribedClass_REF_dec(A);
#ifndef SIMPLE_TEST
DescribedClass_REF_dec(B);
DescribedClass_REF_dec(C);
DescribedClass_REF_dec(D);
#endif /* ! SIMPLE_TEST */

DescribedClass_REF_def(A);
#ifndef SIMPLE_TEST
DescribedClass_REF_def(B);
DescribedClass_REF_def(C);
DescribedClass_REF_def(D);
#endif /* ! SIMPLE_TEST */

main()
{
  ClassDesc::list_all_classes();

  A a;
  cout << "A name:" << a.class_name() << '\n';

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

}
