
#include <comm/picl/picl.h>

ptred2_(a, lda, n, m, p, id, d, e, z, work)
int    *lda,*n,*m,*p,*id;
double *a,*d,*e,*work,*z;
{
  if (*p==1) ptred_single(a, lda, n, m, p, id, d, e, z, work);
  else ptred_parallel(a, lda, n, m, p, id, d, e, z, work);
}

/* ******************************************************** */
/* Function of this subroutine :                            */ 
/*    tridiagonalize a real, symmetric matrix using         */ 
/*    Householder transformation                            */
/* Parameters :                                             */
/*    a[lda][m] - locally held submatrix                    */
/*    lda - leading dimension of array a                    */
/*    n   - size of the matrix a                            */
/*    m   - number of locally held columns                  */
/*    p   - number of nodes used                            */
/*    id  - my node id                                      */
/*  on return :                                             */
/*    d[n]    - the diagonal of the tridiagonal result      */
/*    e[n]    - the offdiagonal of the result(e[1]-e[n-1])  */
/*    z[m][n] - m rows of transformation matrix             */
/*    matrix a will be destroyed                            */
/* -------------------------------------------------------- */

ptred_single(a,lda,n,m,p,id,d,e,z,work)
int    *lda,*n,*m,*p,*id;
double *a,*d,*e,*work,*z;
{
   double  alpha, beta, gamma, dnrm2_(), ddot_(), alpha2; 
   double oobeta,atemp;
   int     i,j,k,l,ld,r,dpsize=sizeof(double);
   int     nproc, slda, sn, sm, sp, sid, q, inc=1;

   /* extract parameters and get  cube information */

   slda = *lda;
   sn = *n;
   sm = *m;
   sp = *p;
   sid = *id;

   /* initialize eigenvector matrix to be identity */

   i = sn * sm;
   alpha2 = 0.0;
   j = 0;
   dcopy_(&i,&alpha2,&j,&z[0],&inc);
   ld = sid;
   for (i=0; i<sm; i++) {
      z[ld*sm+i] = 1.0;
      ld += sp; 
   }

   /* start reduction - one column at a time */

   l = 0;
   ld = sid;
   d[0] = 0.0;
   e[0] = 0.0;
   if (sid == 0) d[0] = a[0];
   for (k=1; k<=sn-1; k++) {
      r = (k-1) % sp;

      /* Use a Householder reflection to zero a(i,k), i = k+2,..., n .*/
      /* Let  a = (0, ..., 0, a(k+1,k) ... a(n,k))',                  */
      /*      u =  a/norm(a) + (k+1-st unit vector),                  */
      /*      beta = -u(k+1) = -norm(u)**2/2 ,                        */
      /*      H = I + u*u'/beta .                                     */
      /* Replace  A  by  H*A*H .                                      */
      /* Store u in D(K+1) through D(N) .                             */
      /* The root node, r, is the owner of column k.                  */

      if (sid == r) {
         q = sn - k;      
         alpha = dnrm2_(&q,&a[l*slda+k],&inc);
         if (a[l*slda+k] < 0.0) alpha = -alpha;
         if (alpha != 0.0) {
            alpha2 = 1.0 / alpha;
            dscal_(&q,&alpha2,&a[l*slda+k],&inc);
            a[l*slda+k] += 1.0;
         }
         dcopy_(&q,&a[l*slda+k],&inc,&d[k],&inc);
         l++;
         ld += sp;
      }

      beta = -d[k];
      if (beta != 0.0) {

         /* Symmetric matrix times vector,  v = A*u.*/
         /* Store  v  in  E(K+1) through E(N) .     */

         alpha2 = 0.0;
         j = 0;
         q = sn - k;
         dcopy_(&q,&alpha2,&j,&e[k],&inc);
         i = ld;
         for (j=l; j<sm; j++) {
            q = sn - i;
            e[i] = e[i] + ddot_(&q,&a[j*slda+i],&inc,&d[i],&inc);
            q = sn - i - 1;
            daxpy_(&q,&d[i],&a[slda*j+i+1],&inc,&e[i+1],&inc);
            i += sp;
         }

         /* v = v/beta            */
         /* gamma = v'*u/(2*beta) */
         /* v = v + gamma*u       */

         if (sid == r) {
            q = sn - k;
            alpha2 = 1.0 / beta;
            dscal_(&q,&alpha2,&e[k],&inc);
            gamma = 0.5*ddot_(&q,&d[k],&inc,&e[k],&inc)/beta;
            daxpy_(&q,&gamma,&d[k],&inc,&e[k],&inc);
         }

         /* Rank two update of A, compute only lower half. */
         /* A  =  A + u'*v + v'*u  =  H*A*H                */

         i = ld;
         for (j=l; j<sm; j++) {
            q = sn - i;
            daxpy_(&q,&d[i],&e[i],&inc,&a[j*slda+i],&inc);
            daxpy_(&q,&e[i],&d[i],&inc,&a[j*slda+i],&inc);
            i += sp;
         }
         q = sn - k;
         oobeta=1.0/beta;
         for (i=0; i<sm; i++) {
            gamma = ddot_(&q,&d[k],&inc,&z[k*sm+i],&sm)*oobeta;
            daxpy_(&q,&gamma,&d[k],&inc,&z[k*sm+i],&sm);
         }
      }

      d[k] = 0.0;
      e[k] = 0.0;
      if (sid == (k % sp)) d[k] = a[l*slda+ld];
      if (sid == r) e[k] = -alpha;
   }
   r = 0;

}

/*
 * Function of this subroutine :
 * tridiagonalize a real, symmetric matrix using
 * Householder transformation
 *
 * Parameters :
 *   a[lda][m] - locally held submatrix
 *   lda - leading dimension of array a
 *   n   - size of the matrix a
 *   m   - number of locally held columns
 *   p   - number of nodes used
 *   id  - my node id
 *
 * on return :
 *   d[n]    - the diagonal of the tridiagonal result
 *   e[n]    - the offdiagonal of the result(e[1]-e[n-1])
 *   z[m][n] - m rows of transformation matrix
 *   matrix a will be destroyed 
 *
 * merge C code from libdmt with FORTRAN code modified by R. Chamberlain
 * FORTRAN COMMENTS:
 *    This version dated 5/4/92
 *    Richard Chamberlain, Intel Supercomputer Systems Division
 *    Improvements:
 *      1. gdcomb of Robert van de Geijn used.
 *      2. look-ahead distribution of Householder vector. Here the node
 *        containing the next Householder vector defers updating the
 *        eigenvector matrix until the next Householder vector is sent.
 */

ptred_parallel(a, lda, n, m, p, id, d, e, z, work)
int *lda, *n, *m, *p, *id;
double *a, *d, *e, *work, *z;
{
  int i, j, k, l, ld, r, dpsize = sizeof(double);
  int kp1l;
  int slda, sn, sm, sp, sid, q, inc = 1;
  double alpha, beta, gamma, dnrm2_(), ddot_(), alpha2;
  double oobeta, atemp;

  /* extract parameters and get  cube information */

  slda = *lda;
  sn = *n;
  sm = *m;
  sp = *p;
  sid = *id;

  /* initialize eigenvector matrix to be identity */

  i = sn * sm;
  alpha2 = 0.0;
  j = 0;
  dcopy_(&i, &alpha2, &j, &z[0], &inc);
  ld = sid;
  for (i = 0; i < sm; i++) {
    z[ld * sm + i] = 1.0;
    ld += sp;
    }

  /* start reduction - one column at a time */

  l = 0;
  ld = sid;
  d[0] = 0.0;
  e[0] = 0.0;
  if (sid == 0) d[0] = a[0];

  for (k = 1; k <= sn - 1; k++) {

    /* Use a Householder reflection to zero a(i,k), i = k+2,..., n .
     * Let  a = (0, ..., 0, a(k+1,k) ... a(n,k))',
     * u =  a/norm(a) + (k+1-st unit vector),
     * beta = -u(k+1) = -norm(u)**2/2,
     * H = I + u*u'/beta.
     * Replace  A  by  H*A*H.
     * Store u in D(K+1) through D(N).
     * The root node, r, is the owner of column k.                  
     */

    r = (k - 1) % sp;
    if (sid == r) {
      kp1l=l*slda+k;
      q = sn - k;
      atemp = a[l * slda + ld];
      alpha = dnrm2_(&q, &a[kp1l], &inc);
      if (a[kp1l] < 0.0) alpha = -alpha;
      if (alpha != 0.0) {
	alpha2 = 1.0 / alpha;
	dscal_(&q, &alpha2, &a[kp1l], &inc);
	a[kp1l] += 1.0;
        }

      bcast0(&a[kp1l], (sn - k) * dpsize, mtype_get() , r);

   /* this is the deferred update of the eigenvector matrix. It was
    * deferred from the last step to accelerate the sending of the Householder
    * vector. Don't do this on the first step.
    */
      if (k != 1) {
	int ik = k - 1; /* ik is a temporary index to the previous step */
	int nmik = sn - ik;

	if (beta != 0.0) {
	  for (i = 0; i < sm; i++) {
	    gamma = ddot_(&nmik, &d[ik], &inc, &z[ik * sm + i], &sm) / beta;
	    daxpy_(&nmik, &gamma, &d[ik], &inc, &z[ik * sm + i], &sm);
	    }
	  }
	e[ik] = 0.0;
	d[ik] = atemp;
        }

   /* now resume normal service */
      dcopy_(&q, &a[kp1l], &inc, &d[k], &inc);
      l++;
      ld += sp;
      }
    else {
      bcast0(&d[k], (sn - k) * dpsize, mtype_get() , r);
      }

    beta = -d[k];
    if (beta != 0.0) {

      /* Symmetric matrix times vector,  v = A*u. */
      /* Store  v  in  E(K+1) through E(N) .     */

      alpha2 = 0.0;
      j = 0;
      q = sn - k;
      dcopy_(&q, &alpha2, &j, &e[k], &inc);
      i = ld;
      for (j = l; j < sm; j++) {
        int ij=j*slda+i;
	q = sn - i;
	e[i] += ddot_(&q, &a[ij], &inc, &d[i], &inc);
	q--;
	daxpy_(&q, &d[i], &a[ij+1], &inc, &e[i + 1], &inc);
	i += sp;
        }
      gop1 (&e[k], sn-k, work, '+', mtype_get());

      /* v = v/beta
       * gamma = v'*u/(2*beta)
       * v = v + gamma*u
       */

      q = sn - k;
      alpha2 = 1.0 / beta;
      dscal_(&q, &alpha2, &e[k], &inc);
      gamma = 0.5 * ddot_(&q, &d[k], &inc, &e[k], &inc) / beta;
      daxpy_(&q, &gamma, &d[k], &inc, &e[k], &inc);

      /* Rank two update of A, compute only lower half. */
      /* A  =  A + u'*v + v'*u  =  H*A*H                */

      i = ld;
      for (j = l; j < sm; j++) {
        double *atmp= &a[j*slda+i];
	q = sn - i;
	daxpy_(&q, &d[i], &e[i], &inc, atmp, &inc);
	daxpy_(&q, &e[i], &d[i], &inc, atmp, &inc);
	i += sp;
        }

      /*  Accumulate m rows of transformation matrix.
       *  Z = Z*H
       *
       * if I have next column, defer updating
       */

      if (sid != k%sp || k == sn - 1) {
	q = sn - k;
	oobeta = 1.0 / beta;
	for (i = 0; i < sm; i++) {
	  gamma = ddot_(&q, &d[k], &inc, &z[k * sm + i], &sm) * oobeta;
	  daxpy_(&q, &gamma, &d[k], &inc, &z[k * sm + i], &sm);
	  }
        }
      }

   /* another bit of calcs to be deferred */
    if (sid != k%sp || k == sn - 1) {
      d[k] = 0.0;
      e[k] = 0.0;
      if (sid == k%sp) d[k] = a[l * slda + ld];
      if (sid == r) e[k] = -alpha;
      }
    }

  /* collect the whole tridiagonal matrix at every node */

  gop1(d, sn, work, '+', mtype_get());
  gop1(e, sn, work, '+', mtype_get());
  }
