
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <util/group/picl.h>

#ifdef HAVE_CUBE_H
#include <cube.h>
#endif
#ifdef HAVE_NX_H
#include <cube.h>
#endif

void
gcollect(double *x, int *lens, double *y)
{
#if defined(HAVE_CUBE_H) || defined(HAVE_NX_H)
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

  leny = 0;
  for (i=0; i < nproc; i++) leny += lens[i];

  memset(y,0,leny);

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
