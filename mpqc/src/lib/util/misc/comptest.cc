
#include <stdio.h>
#include "compute.h"

//Result_dec(int);

class A: public Compute {
  public:
    A();
    Resultint i;
    void compute();
};

A::A():i(this)
{
}

void
A::compute()
{
  i.result_noupdate() = 5;
  i.computed() = 1;
}

main()
{
  A a;

  printf("a.i = %d\n",(int)a.i);
}

  

