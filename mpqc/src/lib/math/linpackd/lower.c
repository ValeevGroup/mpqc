
extern double ddot();
extern double dnrm2();

daxpy_(a,b,c,d,e,f)
void *a;
void *b;
void *c;
void *d;
void *e;
void *f;
{
  daxpy(a,b,c,d,e,f);
  }

dcopy_(a,b,c,d,e)
void *a;
void *b;
void *c;
void *d;
void *e;
{
  dcopy(a,b,c,d,e);
  }

dscal_(a,b,c,d)
void *a;
void *b;
void *c;
void *d;
{
  dscal(a,b,c,d);
  }

dswap_(a,b,c,d,e)
void *a;
void *b;
void *c;
void *d;
void *e;
{
  dswap(a,b,c,d,e);
  }

dsvdc_(a,b,c,d,e,f,g,h,i,j,k,l,m)
void *a;
void *b;
void *c;
void *d;
void *e;
void *f;
void *g;
void *h;
void *i;
void *j;
void *k;
void *l;
void *m;
{
  dsvdc(a,b,c,d,e,f,g,h,i,j,k,l,m);
  }

drotg_(a,b,c,d)
void *a;
void *b;
void *c;
void *d;
{
  drotg(a,b,c,d);
  }

drot_(a,b,c,d,e,f,g)
void *a;
void *b;
void *c;
void *d;
void *e;
void *f;
void *g;
{
  drot(a,b,c,d,e,f,g);
  }

double
ddot_(a,b,c,d,e)
void *a;
void *b;
void *c;
void *d;
void *e;
{
  return ddot(a,b,c,d,e);
  }

double
dnrm2_(a,b,c)
void *a;
void *b;
void *c;
{
  return dnrm2(a,b,c);
  }
