
#include <stdio.h>
#include <util/container/array.h>
#include <util/state/state.h>

main()
{
  int i;
  SSBArray<int> a(10);

  for (i=0; i<10; i++) {
      a(i) = i;
    }

  for (i=0; i<10; i++) {
      printf("a(%d) = %d\n",i,a(i));
    }

  ///////////////
  int j;
  SSBArray2<int> b(3,3);

  for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
          b(i,j) = i + j;
        }
    }

  for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
          printf("b(%d,%d) = %d\n",i,j,b(i,j));
        }
    }
}
