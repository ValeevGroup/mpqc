
 /* The high level PICL routines written in terms of gopf2. */

#include <stdio.h>
#include <stdlib.h>
#include <comm/picl/picl.h>

static int size;

#define gopf_func(opname,op,type) \
  static long opname ## _ ## type (type*x,type*t,int l) \
  { \
  int i; \
  for (i=0; i<size/sizeof(type); i++) x[i] = x[i] op t[i]; \
  return 0; \
  }

#define gopf_xor(type) \
  static long xor_ ## type (type*x,type*t,int l) \
  { \
  int i; \
  for (i=0; i<size/sizeof(type); i++) { \
    if (!(x[i]&&t[i]) && (x[i]||t[i])) x[i] = 1; \
    else x[i] = 0; \
    } \
  return 0; \
  }

#define gopf_minmax(opname,op,type) \
  static long opname ## _ ## type (type*x,type*t,int l) \
  { \
  int i; \
  for (i=0; i<size/sizeof(type); i++) { if (t[i] op x[i]) x[i] = t[i]; }\
  return 0; \
  }

gopf_func(  and, &&, char)
gopf_func(  and, &&, short)
gopf_func(  and, &&, int)
gopf_func(  and, &&, long)
gopf_func(  and, &&, float)
gopf_func(  and, &&, double)
gopf_func(   or, ||, char)
gopf_func(   or, ||, short)
gopf_func(   or, ||, int)
gopf_func(   or, ||, long)
gopf_func(   or, ||, float)
gopf_func(   or, ||, double)
gopf_func( prod,  *, char)
gopf_func( prod,  *, short)
gopf_func( prod,  *, int)
gopf_func( prod,  *, long)
gopf_func( prod,  *, float)
gopf_func( prod,  *, double)
gopf_func(  sum,  +, char)
gopf_func(  sum,  +, short)
gopf_func(  sum,  +, int)
gopf_func(  sum,  +, long)
gopf_func(  sum,  +, float)
gopf_func(  sum,  +, double)
gopf_xor(char)
gopf_xor(short)
gopf_xor(int)
gopf_xor(long)
gopf_xor(float)
gopf_xor(double)
gopf_minmax(  min,  <, char)
gopf_minmax(  min,  <, short)
gopf_minmax(  min,  <, int)
gopf_minmax(  min,  <, long)
gopf_minmax(  min,  <, float)
gopf_minmax(  min,  <, double)
gopf_minmax(  max,  >, char)
gopf_minmax(  max,  >, short)
gopf_minmax(  max,  >, int)
gopf_minmax(  max,  >, long)
gopf_minmax(  max,  >, float)
gopf_minmax(  max,  >, double)

static void do_global(data,len,func,root,msgtype)
void *data;
int len;
long (*func)();
int root;
int msgtype;
{
  char* scratch = (char*) malloc(len);
  size = len;
  gopf2((char*)data,len,scratch,msgtype,root,func);
  free(scratch);
}

void gand0(data, ndata, data_type, msgtype, root)
     void *data;
     int ndata, data_type,msgtype, root;
{
  switch(data_type)
    {
    case 0:
      do_global(data,ndata*sizeof(char),and_char,root,msgtype);
      break;
    case 1:
      do_global(data,ndata*sizeof(short),and_short,root,msgtype);
      break;
    case 2:
      do_global(data,ndata*sizeof(int),and_int,root,msgtype);
      break;
    case 3:
      do_global(data,ndata*sizeof(long),and_long,root,msgtype);
      break;
    case 4:
      do_global(data,ndata*sizeof(float),and_float,root,msgtype);
      break;
    case 5:
      do_global(data,ndata*sizeof(double),and_double,root,msgtype);
      break;
    }
}

void gmax0(data, ndata, data_type, msgtype, root)
     void *data;
     int ndata, data_type,msgtype, root;
{
  switch(data_type)
    {
    case 0:
      do_global(data,ndata*sizeof(char),max_char,root,msgtype);
      break;
    case 1:
      do_global(data,ndata*sizeof(short),max_short,root,msgtype);
      break;
    case 2:
      do_global(data,ndata*sizeof(int),max_int,root,msgtype);
      break;
    case 3:
      do_global(data,ndata*sizeof(long),max_long,root,msgtype);
      break;
    case 4:
      do_global(data,ndata*sizeof(float),max_float,root,msgtype);
      break;
    case 5:
      do_global(data,ndata*sizeof(double),max_double,root,msgtype);
      break;
    }
}

void gmin0(data, ndata, data_type, msgtype, root)
     void *data;
     int ndata, data_type,msgtype, root;
{
  switch(data_type)
    {
    case 0:
      do_global(data,ndata*sizeof(char),min_char,root,msgtype);
      break;
    case 1:
      do_global(data,ndata*sizeof(short),min_short,root,msgtype);
      break;
    case 2:
      do_global(data,ndata*sizeof(int),min_int,root,msgtype);
      break;
    case 3:
      do_global(data,ndata*sizeof(long),min_long,root,msgtype);
      break;
    case 4:
      do_global(data,ndata*sizeof(float),min_float,root,msgtype);
      break;
    case 5:
      do_global(data,ndata*sizeof(double),min_double,root,msgtype);
      break;
    }
}

void gor0(data, ndata, data_type, msgtype, root)
     void *data;
     int ndata, data_type,msgtype, root;
{
  switch(data_type)
    {
    case 0:
      do_global(data,ndata*sizeof(char),or_char,root,msgtype);
      break;
    case 1:
      do_global(data,ndata*sizeof(short),or_short,root,msgtype);
      break;
    case 2:
      do_global(data,ndata*sizeof(int),or_int,root,msgtype);
      break;
    case 3:
      do_global(data,ndata*sizeof(long),or_long,root,msgtype);
      break;
    case 4:
      do_global(data,ndata*sizeof(float),or_float,root,msgtype);
      break;
    case 5:
      do_global(data,ndata*sizeof(double),or_double,root,msgtype);
      break;
    }
}

void gprod0(data, ndata, data_type, msgtype, root)
     void *data;
     int ndata, data_type,msgtype, root;
{
  switch(data_type)
    {
    case 0:
      do_global(data,ndata*sizeof(char),prod_char,root,msgtype);
      break;
    case 1:
      do_global(data,ndata*sizeof(short),prod_short,root,msgtype);
      break;
    case 2:
      do_global(data,ndata*sizeof(int),prod_int,root,msgtype);
      break;
    case 3:
      do_global(data,ndata*sizeof(long),prod_long,root,msgtype);
      break;
    case 4:
      do_global(data,ndata*sizeof(float),prod_float,root,msgtype);
      break;
    case 5:
      do_global(data,ndata*sizeof(double),prod_double,root,msgtype);
      break;
    }
}

void gsum0(data, ndata, data_type, msgtype, root)
     void *data;
     int ndata, data_type,msgtype, root;
{
  switch(data_type)
    {
    case 0:
      do_global(data,ndata*sizeof(char),sum_char,root,msgtype);
      break;
    case 1:
      do_global(data,ndata*sizeof(short),sum_short,root,msgtype);
      break;
    case 2:
      do_global(data,ndata*sizeof(int),sum_int,root,msgtype);
      break;
    case 3:
      do_global(data,ndata*sizeof(long),sum_long,root,msgtype);
      break;
    case 4:
      do_global(data,ndata*sizeof(float),sum_float,root,msgtype);
      break;
    case 5:
      do_global(data,ndata*sizeof(double),sum_double,root,msgtype);
      break;
    }
}

void gxor0(data, ndata, data_type, msgtype, root)
     void *data;
     int ndata, data_type,msgtype, root;
{
  switch(data_type)
    {
    case 0:
      do_global(data,ndata*sizeof(char),xor_char,root,msgtype);
      break;
    case 1:
      do_global(data,ndata*sizeof(short),xor_short,root,msgtype);
      break;
    case 2:
      do_global(data,ndata*sizeof(int),xor_int,root,msgtype);
      break;
    case 3:
      do_global(data,ndata*sizeof(long),xor_long,root,msgtype);
      break;
    case 4:
      do_global(data,ndata*sizeof(float),xor_float,root,msgtype);
      break;
    case 5:
      do_global(data,ndata*sizeof(double),xor_double,root,msgtype);
      break;
    }
}

static void (*comb_routine)(char*p1,char*p2,int,int);
static int comb_routine_datatype;
static int comb_routine_ndata;

static long comb_function(char*x,char*t,int l)
{
  int i;
  comb_routine(x,t,comb_routine_ndata,comb_routine_datatype);
  return 0;
}

void gcomb0(data, ndata, data_type, msgtype, root, comb)
     void *data;
     int ndata, data_type,msgtype, root;
     void (*comb)(char*p1,char*p2,int ndata,int datatype);
{
  comb_routine = comb;
  comb_routine_datatype = data_type;
  comb_routine_ndata = ndata;
  switch(data_type)
    {
    case 0:
      do_global(data,ndata*sizeof(char),comb_function,root,msgtype);
      break;
    case 1:
      do_global(data,ndata*sizeof(short),comb_function,root,msgtype);
      break;
    case 2:
      do_global(data,ndata*sizeof(int),comb_function,root,msgtype);
      break;
    case 3:
      do_global(data,ndata*sizeof(long),comb_function,root,msgtype);
      break;
    case 4:
      do_global(data,ndata*sizeof(float),comb_function,root,msgtype);
      break;
    case 5:
      do_global(data,ndata*sizeof(double),comb_function,root,msgtype);
      break;
    }
}

