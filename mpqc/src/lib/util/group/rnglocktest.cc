
#define TEST 1
#define VERBOSE 0

#include <math.h>
#include <util/group/rnglock.h>

main()
{
  const int size = 10000;
  int bin[size];
  int i;

  RangeLock lock;

  for (i=0; i<size; i++) bin[i] = 0;

  for (i=0; i<1000000; i++) {
      int start = random()%size;
      int length = random()%30;
      if (length == 0) length++;
      int fence = start + length;
      if (fence > size) fence = size;
      int val = random()%2 ? -1:1;
      val = 1;
#if VERBOSE
      printf("adding block %d: start = %d, fence = %d, val = %d\n",
             i, start, fence, val);
#endif
      lock.sum(start, fence, val);
      for (int j=start; j<fence; j++) {
          bin[j] += val;
        }
    }

  printf("printing nonzero deltas\n");

  for (i=0; i<size; i++) {
      int delta = bin[i] - lock.lockvalue(i);
      if (delta) printf(" %5d: %8d (%8d %8d)\n", i, delta,
                        bin[i], lock.lockvalue(i));
    }

  return 0;
}
