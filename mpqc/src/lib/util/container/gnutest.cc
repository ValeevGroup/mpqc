
#include <util/container/set.h>

int
main()
{
  int i;
  AVLSet<int> iset;

  for (i=10; i<20; i++) {
      printf("adding %d\n",i);
      iset.add(i);
    }

  for (i=1; i<15; i++) {
      printf("adding %d\n",i);
      iset.add(i);
    }

  int ten=10, fifteen=15;
  iset.del(ten);
  iset.del(fifteen);

  printf("iset:");
  Pix I;
  for (I=iset.first(); I; iset.next(I)) {
      printf(" %d",iset(I));
    }
  printf("\n");

  ///////////////////////////////////////////////////////////////

  Arrayset<int> aset;

  for (I=aset.first(); I; aset.next(I)) {
      printf(" %d",aset(I));
    }

  for (i=10; i<20; i++) {
      printf("adding %d\n",i);
      aset.add(i);
    }

  for (i=1; i<15; i++) {
      printf("adding %d\n",i);
      aset.add(i);
    }

  aset.del(ten);
  aset.del(fifteen);

  printf("aset:");
  for (I=aset.first(); I; aset.next(I)) {
      printf(" %d",aset(I));
    }
  printf("\n");

  return 0;
}
