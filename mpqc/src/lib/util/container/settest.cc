
#include <stdio.h>
#include "ref.h"
#include "array.h"
#include "set.h"

class Intg: public RefCount {
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
void Intg::print() { printf("%d\n",_int); }

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
  RefB b;
  //RefA a(b); // illegal
  RefA a((A*)b);
  
  Intg* i1 = new Intg(101);
  i1->print();

  RefIntg ii = new Intg(100);

  ArraysetRefIntg as;

  for (int i=0; i<2; i++) {
      as.add(new Intg(i));
    }
  as[1] = ii;
  for (i=0; i<2; i++) {
      as[i]->print();
    }
  ii->print();
  i1->print();
}
