
#include <stdio.h>
#include <util/misc/compute.h>

#ifdef __GNUC__
template class NCResult<int>;
#endif

class A: public Compute {
  public:
    A();
    Resultint i;
    AccResultdouble a;
    void compute();
    void print();
};

A::A():i(this),a(this)
{
  a.set_desired_accuracy(0.1);
  a.set_actual_accuracy(0.1);
  a.result_noupdate() = 0.0;
  a.computed() = 1;
}

void
A::compute()
{
  printf("computing");
  if (i.needed()) {
      i.result_noupdate() = 5;
      i.computed() = 1;
      printf(" i");
    }
  if (a.needed()) {
      a.result_noupdate() += 0.001;
      a.computed() = 1;
      printf(" a");
      a.set_actual_accuracy(a.desired_accuracy());
    }
  printf("\n");
}

void
A::print()
{
  printf("A: i = %d, a = %5.3f\n",(int)i,(double)a);
}

main()
{
  A a;

  printf("should not compute a\n");
  a.print();

  a.a.set_desired_accuracy(0.01);

  printf("should compute a\n");
  a.print();

  a.a.set_desired_accuracy(0.1);

  printf("should not compute a\n");
  a.print();

  a.a.set_desired_accuracy(0.01);

  printf("should not compute a\n");
  a.print();

  a.a.set_desired_accuracy(0.001);

  printf("should compute a\n");
  a.print();
}
