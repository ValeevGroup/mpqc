
extern double DDOT();
extern double DNRM2();

daxpy_(a,b,c,d,e,f)
void *a;
void *b;
void *c;
void *d;
void *e;
void *f;
{
  DAXPY(a,b,c,d,e,f);
  }

dcopy_(a,b,c,d,e)
void *a;
void *b;
void *c;
void *d;
void *e;
{
  DCOPY(a,b,c,d,e);
  }

dscal_(a,b,c,d)
void *a;
void *b;
void *c;
void *d;
{
  DSCAL(a,b,c,d);
  }

dswap_(a,b,c,d,e)
void *a;
void *b;
void *c;
void *d;
void *e;
{
  DSWAP(a,b,c,d,e);
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
  DSVDC(a,b,c,d,e,f,g,h,i,j,k,l,m);
  }

drotg_(a,b,c,d)
void *a;
void *b;
void *c;
void *d;
{
  DROTG(a,b,c,d);
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
  DROT(a,b,c,d,e,f,g);
  }

double
ddot_(a,b,c,d,e)
void *a;
void *b;
void *c;
void *d;
void *e;
{
  return DDOT(a,b,c,d,e);
  }

double
dnrm2_(a,b,c)
void *a;
void *b;
void *c;
{
  return DNRM2(a,b,c);
  }
