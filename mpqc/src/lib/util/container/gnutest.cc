
#include "set.h"

int
main()
{
  AVLSet<int> iset;

  for (int i=10; i<20; i++) {
      printf("adding %d\n",i);
      iset.add(i);
    }

  for (i=1; i<15; i++) {
      printf("adding %d\n",i);
      iset.add(i);
    }

  iset.del(10);
  iset.del(15);

  printf("iset:");
  for (Pix I=iset.first(); I; iset.next(I)) {
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

  aset.del(10);
  aset.del(15);

  printf("aset:");
  for (I=aset.first(); I; aset.next(I)) {
      printf(" %d",aset(I));
    }
  printf("\n");

  return 0;
}
