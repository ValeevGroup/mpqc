
#include <stdio.h>
#include "ref.h"
#include "array.h"
#include "set.h"

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
void Intg::print() { printf("%d\n",_int); }

REF_dec(Intg);
REF_def(Intg);

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
  printf("0x%08x ref count = %d\n",aap,aap->nreference());
//   delete aap;
//   printf("0x%08x ref count = %d\n",aap,aap->_reference_count_);
  aar = aap;
//   delete aap;
//   printf("0x%08x ref count = %d\n",aap,aap->_reference_count_);
  aar = 0;
  
  RefB b;
  //RefA a(b); // illegal
  RefA a(b.pointer());
  
  Intg* i1 = new Intg(101);
  i1->print();

  RefIntg ii = new Intg(100);

  Arrayset<RefIntg> as;

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
