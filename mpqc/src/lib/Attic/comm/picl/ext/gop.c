#include <stdio.h>
#include <stdlib.h>
#include <comm/picl/picl.h>
#include "piclext.h"

#if defined(PARAGON)
#include <nx.h>
#endif

#include "assert.h"

void gop1(val, len, tmp, op, type)
double *val;
int op;
double *tmp;
int len, type;
{

#if defined(PARAGON)
  switch (op) {
    case 'm':
      gdlow(val,len,tmp);
      break;
    case 'M':
      gdhigh(val,len,tmp);
      break;
    case '+':
      gdsum(val,len,tmp);
      break;
    }
#elif defined(I860)
  switch (op) {
    case 'm':
      gmin0(val,len,5,type,0);
      bcast0(val,len,type,0);
      break;
    case 'M':
      gmax0(val,len,5,type,0);
      bcast0(val,len,type,0);
      break;
    case '+':
      gdcomb(len,val,tmp,0,cubedim0());
      break;
    }
#else
  switch (op) {
    case 'm':
      gmin0(val,len,5,type,0);
      break;
    case 'M':
      gmax0(val,len,5,type,0);
      break;
    case '+':
      gsum0(val,len,5,type,0);
      break;
    }
  bcast0(val,len*sizeof(double),type,0);
#endif

}

void gop0(val,len,op,type)
double *val;
int op;
int len, type;
{
  double *tmp = (double*)malloc(sizeof(double)*len);
  if (!tmp) {
      fprintf(stderr,"gop0: couldn't allocate %d bytes\n",len);
      exit(1);
    }
  gop1(val,len,tmp,op,type);
  free(tmp);
}

/* These routines are used to implement a gop0 for signed chars. */
static void min_schar (x,t,l,type)
signed char*x;signed char*t;int l;int type;
{
  int i;
  for (i=0; i<l; i++) { if (t[i] < x[i]) x[i] = t[i]; }
}
static void max_schar (x,t,l,type)
signed char*x;signed char*t;int l;int type;
{
  int i;
  for (i=0; i<l; i++) { if (t[i] > x[i]) x[i] = t[i]; }
}
static void sum_schar (x,t,l,type)
signed char*x;signed char*t;int l;int type;
{
  int i;
  for (i=0; i<l; i++) x[i] += t[i];
}

/* A gop0 for signed chars. */
void gop0_sc(signed char*val,int len,int op,int type)
{
  switch (op) {
    case 'm':
      gcomb0(val,len,0,type,0,min_schar);
      break;
    case 'M':
      gcomb0(val,len,0,type,0,max_schar);
      break;
    case '+':
      gcomb0(val,len,0,type,0,sum_schar);
      break;
    }
  bcast0(val,len,type,0);
}

