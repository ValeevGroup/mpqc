
#include <stdio.h>
#include <util/container/ref.h>
#include <util/container/reftestx.h>

int
main()
{
  {
  Ref<X> x1 = new X;
  RefX x2 = new X;
#if REF_MANAGE
  x2->unmanage();
#endif
#if REF_CHECK_STACK
    {
      X z;
      RefX zz(&z);
    }
#endif

  x2 = x1.pointer();
  if (x1 != x1) abort();
  if (x2 != x2) abort();

  for (int i=1000000; i; i--) {
      x1->reference();
    }
  for (i=1000000; i; i--) {
      x1->dereference();
    }
  for (i=1000000; i; i--) {
      Ref<X> y = x1;
    }

  fprintf(stderr,"nx = %d (inner scope)\n",X::nx);
  }

  fprintf(stderr,"nx = %d (outer scope)\n",X::nx);
  return 0;
}

