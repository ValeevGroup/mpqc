
#include <stdio.h>
#include <util/group/picl.h>

int
prober(start,end)
int start;
int end;
{
  int i;
  int found=0;

  sync0();

  for (i=start; i<=end; i++) {
    if (probe0(i)) {
      found = 1;
      printf("found message of type %d on node %d\n",i,mynode0());
      }
    }

  sync0();
  return found;
  }

