
#include <stdio.h>
#include <stdlib.h>

#include <util/group/picl.h>

#if defined(I860) && !defined(PARAGON)
#include <cube.h>
#endif

void
gcollect(double *x, int *lens, double *y)
{
#if defined(I860) && !defined(PARAGON)
  gcolx(x,lens,y);
#else

  int i,leny;
  int me = mynode0();
  int nproc = numnodes0();
  int dim = cubedim0();
  int myoff=0;
  double *temp;

  if (dim==0) {
    memcpy(y,x,lens[0]);
    return;
    }
    
  for (i=0; i < me; i++) myoff += lens[i];
  myoff /= sizeof(double);

  for (i=0; i < nproc; i++) leny += lens[i];
  bzero(y,leny);

  memcpy(&y[myoff],x,lens[me]);

  temp = (double *) malloc(leny);
  if (!temp) {
    fprintf(stderr,"gcollect: could not malloc temp %d\n",leny);
    exit(1);
    }

  gop1(y,leny/sizeof(double),temp,'+',3);

  free(temp);
#endif
  }
