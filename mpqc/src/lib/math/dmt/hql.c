/* ******************************************************** */
/* Function of this subroutine :                            */ 
/*    Diagonalize a real, symmetric matrix                  */ 
/*                                                          */
/* Parameters :                                             */
/*    tci - channel number (not used)                       */
/*    n   - size of the matrix                              */
/*    m   - number of locally held columns                  */
/*    a[n][m] - locally held submatrix                      */
/*    d[n]    - returned with eigenvalues                   */
/*    e[n]    - scratch space                               */
/*    sigma[n]- scratch space                               */
/*    z[m][n] - scratch space                               */
/*    v[n][m] - returned with eigenvectors                  */
/*    w[3*n]  - scratch space                               */
/*    ind[n]  - scratch space (integer)                     */
/* -------------------------------------------------------- */

#include <stdio.h>
#include <math.h>

#include <tmpl.h>
#include <util/misc/libmisc.h>
#include <util/group/picl.h>

static double absol();
static double epslon ();
static void update();
static void pimtql2_ ();
static void pflip();

extern void ptred2_();
extern void dcopy_();
extern void dswap_();

 
void
diagonal_(tci,tn,tm,a,d,e,sigma,z,v,w,ind)
int      *tci, *tn, *tm, *ind;
double   *a, *d, *e, *sigma, *z, *v, *w;
{
  int m,n,id,i,info,one=1;
  int nproc,host;

  who0(&nproc,&id,&host);

  n = *tn;
  m = *tm;

 /* reduce A to tridiagonal matrix using Householder transformation */

  ptred2_(&a[0],&n,&n,&m,&nproc,&id,&d[0],&e[0],&z[0],&w[0]);

 /* diagonalize tridiagonal matrix using implicit QL method */

  pimtql2_(d,e,&n,z,&m,&info);
  if (info != 0) message0("Nonzero ierr from psytqr");

 /* rearrange the eigenvectors by transposition */

  i = m * n;
  dcopy_(&i,&z[0],&one,&a[0],&one);
  pflip(id,n,m,nproc,&a[0],&v[0],&w[0]);
}

/* ******************************************************** */
/* Function of this subroutine :                            */ 
/*    diagonalize a real, symmetric tridiagonal matrix      */
/*    using the QL method                                   */ 
/* Parameters :                                             */
/*  on entry :                                              */
/*    d[n]    - the diagonal of the tridiagonal result      */
/*    e[n]    - the offdiagonal of the result(e[1]-e[n-1])  */
/*    sn      - size of the tridiagonal matrix              */
/*    z[m][n] - m rows of transformation matrix from before */
/*    m   - number of locally held columns                  */
/*  on return :                                             */
/*    d[n]    - the eigenvalues                             */
/*    e[n]    - non-useful information                      */
/*    z[m][n] - m rows of eigenvectors                      */
/*    info    - if 0, results are accurate                  */
/*              if nonzero, results may be inaccurate       */
/* -------------------------------------------------------- */

static void
pimtql2_ (d,e,sn,z,sm,info)
int    *sm,*sn,*info;
double *d,*e,*z;
{
   double  c,s,t,q,u,p,h,macheps;
   int     n,m,i,j,k,im,its,maxit=30,one=1;

   /* extract parameters */

   *info = 0;
   n = *sn;
   m = *sm;
   t = 1.0;
   macheps = epslon(t);

   for (i=1; i<n; i++) e[i-1] = e[i];
   e[n-1] = 0.0;
   k = n - 2;
   for (j=0; j<n; j++) {
      its = 0;
      while (its < maxit) {
         for (im=j; im<=k; im++) 
            if (absol(e[im]) <= macheps*(absol(d[im])+absol(d[im+1]))) break; 
         u = d[j];
         if (im == j) break;
         else {
            its++;

            /* form implicit Wilkinson shift */

            q = (d[j+1] - u) / (2.0 * e[j]);
            t = sqrt(1.0 + q * q);      
            q = d[im] - u + e[j]/((q < 0.0) ? q - t : q + t);
            u = 0.0;
            s = c = 1.0;
            for (i=im-1; i>=j; i--) {
               p = s * e[i];
               h = c * e[i];
               if (absol(p) >= absol(q)) {
                  c = q / p;
                  t = sqrt(1.0 + c * c);
                  e[i+1] = p * t;
                  s = 1.0 / t;
                  c *= s;
               } else {
                  s = p / q;
                  t = sqrt(1.0 + s * s);
                  e[i+1] = q * t;
                  c = 1.0 / t;
                  s *= c;
               }
               q = d[i+1] - u;
               t = (d[i] - q) * s + 2.0 * c * h;
               u = s * t;
               d[i+1] = q + u;
               q = c * t - h;

               /* form eigenvectors */

#if 0
               for (ia=0; ia<m; ia++) {
                  p = z[(i+1)*m+ia];
                  z[(i+1)*m+ia] = s * z[i*m+ia] + c * p; 
                  z[i*m+ia] = c * z[i*m+ia] - s * p;
               }
#else
               update(&z[i*m],m,c,s);
#endif
            }
            d[j] -= u;
            e[j] = q;
            e[im] = 0.0;
         }
      }
      if (its == maxit) {
         *info = its;
         break;
      }
   }

   /* order eigenvalues and eigenvectors */

   for (j=0; j<n-1; j++) {
      k = j;
      for (i=j+1; i<n; i++) if (d[i] < d[k]) k = i;
      if (k != j) {
         dswap_(&one,&d[j],&one,&d[k],&one);
         dswap_(&m,&z[j*m],&one,&z[k*m],&one); 
      }
   }
}

/* ******************************************************** */

static double 
absol(x)
double x;
{
 if (x > 0.0)
   return(x);
 else
   return(-x);
}

/* ******************************************************** */
/* Function : transpose matrix                              */
/* -------------------------------------------------------- */

static void
pflip(id,n,m,p,ar,ac,w)
int    id,n,m,p;
double *ar,*ac,*w;
{
  int i,k,r,dpsize=sizeof(double),one=1;

  i = 0;
  for (k=0; k<n; k++) {
    r = k % p;
    if (id == r) {
      dcopy_(&n,&ar[i],&m,&w[0],&one);
      i++;
    }
    bcast0 (&w[0], n*dpsize, mtype_get(), r);
    dcopy_(&m,&w[id],&p,&ac[k],&n);
  }
}

/* ******************************************************** */
/* Function : calculate machine epslon                      */
/* -------------------------------------------------------- */

static double
epslon (x) 
double x; 
{
  double a,b,c,eps; 

  a = 4.0/3.0;
  eps = 0.0;
  while (eps == 0.0) {
    b = a - 1.0; 
    c = 3.0 * b; 
    eps = c-1.0; 
    if (eps < 0.0) eps = - eps;
  }
  if (x < 0.0) a = - x;
  return(eps*a); 
}

static void
update(z,m,c,s)
double *z;
int m;
double c;
double s;
{
  register int i;
  register double p;

  for (i=0; i < m; i++) {
    p = z[i+m];
    z[m+i] = s * z[i] + c * p;
    z[i]   = c * z[i] - s * p;
  }
}
