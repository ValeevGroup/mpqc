
/* Modified by c janssen to be a standalone C program. */

/* linpackd.f -- translated by f2c (version 19940305).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <math.h>
typedef double doublereal;
typedef double real;
typedef int integer;
typedef int logical;
#define TRUE_ 1
#define FALSE_ 0
#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

static double
d_sign(a,b)
    doublereal *a, *b;
{
  double x;
  x = (*a >= 0 ? *a : - *a);
  return( *b >= 0 ? x : -x);
}

/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b90 = -1.;
static doublereal c_b149 = 1.;

/* Subroutine */ int daxpy_(n, da, dx, incx, dy, incy)
integer *n;
doublereal *da, *dx;
integer *incx;
doublereal *dy;
integer *incy;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m, ix, iy, mp1;


/*     CONSTANT TIMES A VECTOR PLUS A VECTOR. */
/*     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*da == 0.) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS */
/*          NOT EQUAL TO 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	dy[iy] += *da * dx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        CODE FOR BOTH INCREMENTS EQUAL TO 1 */


/*        CLEAN-UP LOOP */

L20:
    m = *n % 4;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= i__1; ++i) {
	dy[i] += *da * dx[i];
/* L30: */
    }
    if (*n < 4) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 4) {
	dy[i] += *da * dx[i];
	dy[i + 1] += *da * dx[i + 1];
	dy[i + 2] += *da * dx[i + 2];
	dy[i + 3] += *da * dx[i + 3];
/* L50: */
    }
    return 0;
} /* daxpy_ */

/* Subroutine */ int dcopy_(n, dx, incx, dy, incy)
integer *n;
doublereal *dx;
integer *incx;
doublereal *dy;
integer *incy;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m, ix, iy, mp1;


/*     COPIES A VECTOR, X, TO A VECTOR, Y. */
/*     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS */
/*          NOT EQUAL TO 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	dy[iy] = dx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        CODE FOR BOTH INCREMENTS EQUAL TO 1 */


/*        CLEAN-UP LOOP */

L20:
    m = *n % 7;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= i__1; ++i) {
	dy[i] = dx[i];
/* L30: */
    }
    if (*n < 7) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 7) {
	dy[i] = dx[i];
	dy[i + 1] = dx[i + 1];
	dy[i + 2] = dx[i + 2];
	dy[i + 3] = dx[i + 3];
	dy[i + 4] = dx[i + 4];
	dy[i + 5] = dx[i + 5];
	dy[i + 6] = dx[i + 6];
/* L50: */
    }
    return 0;
} /* dcopy_ */

doublereal ddot_(n, dx, incx, dy, incy)
integer *n;
doublereal *dx;
integer *incx;
doublereal *dy;
integer *incy;
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i, m;
    static doublereal dtemp;
    static integer ix, iy, mp1;


/*     FORMS THE DOT PRODUCT OF TWO VECTORS. */
/*     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS */
/*          NOT EQUAL TO 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	dtemp += dx[ix] * dy[iy];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    ret_val = dtemp;
    return ret_val;

/*        CODE FOR BOTH INCREMENTS EQUAL TO 1 */


/*        CLEAN-UP LOOP */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= i__1; ++i) {
	dtemp += dx[i] * dy[i];
/* L30: */
    }
    if (*n < 5) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 5) {
	dtemp = dtemp + dx[i] * dy[i] + dx[i + 1] * dy[i + 1] + dx[i + 2] * 
		dy[i + 2] + dx[i + 3] * dy[i + 3] + dx[i + 4] * dy[i + 4];
/* L50: */
    }
L60:
    ret_val = dtemp;
    return ret_val;
} /* ddot_ */

doublereal dnrm2_(n, dx, incx)
integer *n;
doublereal *dx;
integer *incx;
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal one = 1.;
    static doublereal cutlo = 8.232e-11;
    static doublereal cuthi = 1.304e19;

    /* Format strings */
    static char fmt_30[] = "";
    static char fmt_50[] = "";
    static char fmt_70[] = "";
    static char fmt_110[] = "";

    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1;

    /* Local variables */
    static doublereal xmax;
    static integer next, i, j, nn;
    static doublereal hitest, sum;

    /* Assigned format variables */
    char *next_fmt;

    /* Parameter adjustments */
    --dx;

    /* Function Body */

/*     EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX() WITH STORAGE */
/*     INCREMENT INCX . */
/*     IF    N .LE. 0 RETURN WITH RESULT = 0. */
/*     IF N .GE. 1 THEN INCX MUST BE .GE. 1 */

/*           C.L.LAWSON, 1978 JAN 08 */

/*     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE */
/*     HOPEFULLY APPLICABLE TO ALL MACHINES. */
/*         CUTLO = MAXIMUM OF  DSQRT(U/EPS)  OVER ALL KNOWN MACHINES. */
/*         CUTHI = MINIMUM OF  DSQRT(V)      OVER ALL KNOWN MACHINES. */
/*     WHERE */
/*         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1. */
/*         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT) */
/*         V   = LARGEST  NO.            (OVERFLOW  LIMIT) */

/*     BRIEF OUTLINE OF ALGORITHM.. */

/*     PHASE 1    SCANS ZERO COMPONENTS. */
/*     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO */
/*     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO */
/*     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M */
/*     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX. */

/*     VALUES FOR CUTLO AND CUTHI.. */
/*     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER */
/*     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS.. */
/*     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE 
*/
/*                   UNIVAC AND DEC AT 2**(-103) */
/*                   THUS CUTLO = 2**(-51) = 4.44089E-16 */
/*     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC. */
/*                   THUS CUTHI = 2**(63.5) = 1.30438E19 */
/*     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC. */
/*                   THUS CUTLO = 2**(-33.5) = 8.23181D-11 */
/*     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19 */
/*     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 / */
/*     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 / */

    if (*n > 0) {
	goto L10;
    }
    ret_val = zero;
    goto L300;

L10:
    next = 0;
    next_fmt = fmt_30;
    sum = zero;
    nn = *n * *incx;
/*                                                 BEGIN MAIN LOOP */
    i = 1;
L20:
    switch ((int)next) {
	case 0: goto L30;
	case 1: goto L50;
	case 2: goto L70;
	case 3: goto L110;
    }
L30:
    if ((d__1 = dx[i], fabs(d__1)) > cutlo) {
	goto L85;
    }
    next = 1;
    next_fmt = fmt_50;
    xmax = zero;

/*                        PHASE 1.  SUM IS ZERO */

L50:
    if (dx[i] == zero) {
	goto L200;
    }
    if ((d__1 = dx[i], fabs(d__1)) > cutlo) {
	goto L85;
    }

/*                                PREPARE FOR PHASE 2. */
    next = 2;
    next_fmt = fmt_70;
    goto L105;

/*                                PREPARE FOR PHASE 4. */

L100:
    i = j;
    next = 3;
    next_fmt = fmt_110;
    sum = sum / dx[i] / dx[i];
L105:
    xmax = (d__1 = dx[i], fabs(d__1));
    goto L115;

/*                   PHASE 2.  SUM IS SMALL. */
/*                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW. */

L70:
    if ((d__1 = dx[i], fabs(d__1)) > cutlo) {
	goto L75;
    }

/*                     COMMON CODE FOR PHASES 2 AND 4. */
/*                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW. 
*/

L110:
    if ((d__1 = dx[i], fabs(d__1)) <= xmax) {
	goto L115;
    }
/* Computing 2nd power */
    d__1 = xmax / dx[i];
    sum = one + sum * (d__1 * d__1);
    xmax = (d__1 = dx[i], fabs(d__1));
    goto L200;

L115:
/* Computing 2nd power */
    d__1 = dx[i] / xmax;
    sum += d__1 * d__1;
    goto L200;


/*                  PREPARE FOR PHASE 3. */

L75:
    sum = sum * xmax * xmax;


/*     FOR REAL OR D.P. SET HITEST = CUTHI/N */
/*     FOR COMPLEX      SET HITEST = CUTHI/(2*N) */

L85:
    hitest = cuthi / (real) (*n);

/*                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING. */

    i__1 = nn;
    i__2 = *incx;
    for (j = i; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
	if ((d__1 = dx[j], fabs(d__1)) >= hitest) {
	    goto L100;
	}
/* L95: */
/* Computing 2nd power */
	d__1 = dx[j];
	sum += d__1 * d__1;
    }
    ret_val = sqrt(sum);
    goto L300;

L200:
    i += *incx;
    if (i <= nn) {
	goto L20;
    }

/*              END OF MAIN LOOP. */

/*              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING. */

    ret_val = xmax * sqrt(sum);
L300:
    return ret_val;
} /* dnrm2_ */

/* Subroutine */ int dscal_(n, da, dx, incx)
integer *n;
doublereal *da, *dx;
integer *incx;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i, m, nincx, mp1;


/*     SCALES A VECTOR BY A CONSTANT. */
/*     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        CODE FOR INCREMENT NOT EQUAL TO 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i = 1; i__2 < 0 ? i >= i__1 : i <= i__1; i += i__2) {
	dx[i] = *da * dx[i];
/* L10: */
    }
    return 0;

/*        CODE FOR INCREMENT EQUAL TO 1 */


/*        CLEAN-UP LOOP */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__2 = m;
    for (i = 1; i <= i__2; ++i) {
	dx[i] = *da * dx[i];
/* L30: */
    }
    if (*n < 5) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i = mp1; i <= i__2; i += 5) {
	dx[i] = *da * dx[i];
	dx[i + 1] = *da * dx[i + 1];
	dx[i + 2] = *da * dx[i + 2];
	dx[i + 3] = *da * dx[i + 3];
	dx[i + 4] = *da * dx[i + 4];
/* L50: */
    }
    return 0;
} /* dscal_ */

/* Subroutine */ int dswap_(n, dx, incx, dy, incy)
integer *n;
doublereal *dx;
integer *incx;
doublereal *dy;
integer *incy;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m;
    static doublereal dtemp;
    static integer ix, iy, mp1;


/*     INTERCHANGES TWO VECTORS. */
/*     USES UNROLLED LOOPS FOR INCREMENTS EQUAL ONE. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL */
/*         TO 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	dtemp = dx[ix];
	dx[ix] = dy[iy];
	dy[iy] = dtemp;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*       CODE FOR BOTH INCREMENTS EQUAL TO 1 */


/*       CLEAN-UP LOOP */

L20:
    m = *n % 3;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= i__1; ++i) {
	dtemp = dx[i];
	dx[i] = dy[i];
	dy[i] = dtemp;
/* L30: */
    }
    if (*n < 3) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 3) {
	dtemp = dx[i];
	dx[i] = dy[i];
	dy[i] = dtemp;
	dtemp = dx[i + 1];
	dx[i + 1] = dy[i + 1];
	dy[i + 1] = dtemp;
	dtemp = dx[i + 2];
	dx[i + 2] = dy[i + 2];
	dy[i + 2] = dtemp;
/* L50: */
    }
    return 0;
} /* dswap_ */

/* Subroutine */ int dsvdc_(x, ldx, n, p, s, e, u, ldu, v, ldv, work, job, 
	info)
doublereal *x;
integer *ldx, *n, *p;
doublereal *s, *e, *u;
integer *ldu;
doublereal *v;
integer *ldv;
doublereal *work;
integer *job, *info;
{
    /* System generated locals */
    integer x_dim1, x_offset, u_dim1, u_offset, v_dim1, v_offset, i__1, i__2, 
	    i__3;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7;

    /* Local variables */
    static integer kase;
    extern doublereal ddot_();
    static integer jobu, iter;
    extern /* Subroutine */ int drot_();
    static doublereal test;
    extern doublereal dnrm2_();
    static integer nctp1;
    static doublereal b, c;
    static integer nrtp1;
    static doublereal f, g;
    static integer i, j, k, l, m;
    static doublereal t, scale;
    extern /* Subroutine */ int dscal_();
    static doublereal shift;
    extern /* Subroutine */ int dswap_(), drotg_();
    static integer maxit;
    extern /* Subroutine */ int daxpy_();
    static logical wantu, wantv;
    static doublereal t1, ztest, el;
    static integer kk;
    static doublereal cs;
    static integer ll, mm, ls;
    static doublereal sl;
    static integer lu;
    static doublereal sm, sn;
    static integer lm1, mm1, lp1, mp1, nct, ncu, lls, nrt;
    static doublereal emm1, smm1;



/*     DSVDC IS A SUBROUTINE TO REDUCE A DOUBLE PRECISION NXP MATRIX X */
/*     BY ORTHOGONAL TRANSFORMATIONS U AND V TO DIAGONAL FORM.  THE */
/*     DIAGONAL ELEMENTS S(I) ARE THE SINGULAR VALUES OF X.  THE */
/*     COLUMNS OF U ARE THE CORRESPONDING LEFT SINGULAR VECTORS, */
/*     AND THE COLUMNS OF V THE RIGHT SINGULAR VECTORS. */

/*     ON ENTRY */

/*         X         DOUBLE PRECISION(LDX,P), WHERE LDX.GE.N. */
/*                   X CONTAINS THE MATRIX WHOSE SINGULAR VALUE */
/*                   DECOMPOSITION IS TO BE COMPUTED.  X IS */
/*                   DESTROYED BY DSVDC. */

/*         LDX       INTEGER. */
/*                   LDX IS THE LEADING DIMENSION OF THE ARRAY X. */

/*         N         INTEGER. */
/*                   N IS THE NUMBER OF COLUMNS OF THE MATRIX X. */

/*         P         INTEGER. */
/*                   P IS THE NUMBER OF ROWS OF THE MATRIX X. */

/*         LDU       INTEGER. */
/*                   LDU IS THE LEADING DIMENSION OF THE ARRAY U. */
/*                   (SEE BELOW). */

/*         LDV       INTEGER. */
/*                   LDV IS THE LEADING DIMENSION OF THE ARRAY V. */
/*                   (SEE BELOW). */

/*         WORK      DOUBLE PRECISION(N). */
/*                   WORK IS A SCRATCH ARRAY. */

/*         JOB       INTEGER. */
/*                   JOB CONTROLS THE COMPUTATION OF THE SINGULAR */
/*                   VECTORS.  IT HAS THE DECIMAL EXPANSION AB */
/*                   WITH THE FOLLOWING MEANING */

/*                        A.EQ.0    DO NOT COMPUTE THE LEFT SINGULAR */
/*                                  VECTORS. */
/*                        A.EQ.1    RETURN THE N LEFT SINGULAR VECTORS */
/*                                  IN U. */
/*                        A.GE.2    RETURN THE FIRST MIN(N,P) SINGULAR */
/*                                  VECTORS IN U. */
/*                        B.EQ.0    DO NOT COMPUTE THE RIGHT SINGULAR */
/*                                  VECTORS. */
/*                        B.EQ.1    RETURN THE RIGHT SINGULAR VECTORS */
/*                                  IN V. */

/*     ON RETURN */

/*         S         DOUBLE PRECISION(MM), WHERE MM=MIN(N+1,P). */
/*                   THE FIRST MIN(N,P) ENTRIES OF S CONTAIN THE */
/*                   SINGULAR VALUES OF X ARRANGED IN DESCENDING */
/*                   ORDER OF MAGNITUDE. */

/*         E         DOUBLE PRECISION(P). */
/*                   E ORDINARILY CONTAINS ZEROS.  HOWEVER SEE THE */
/*                   DISCUSSION OF INFO FOR EXCEPTIONS. */

/*         U         DOUBLE PRECISION(LDU,K), WHERE LDU.GE.N.  IF */
/*                                   JOBA.EQ.1 THEN K.EQ.N, IF JOBA.GE.2 
*/
/*                                   THEN K.EQ.MIN(N,P). */
/*                   U CONTAINS THE MATRIX OF RIGHT SINGULAR VECTORS. */
/*                   U IS NOT REFERENCED IF JOBA.EQ.0.  IF N.LE.P */
/*                   OR IF JOBA.EQ.2, THEN U MAY BE IDENTIFIED WITH X */
/*                   IN THE SUBROUTINE CALL. */

/*         V         DOUBLE PRECISION(LDV,P), WHERE LDV.GE.P. */
/*                   V CONTAINS THE MATRIX OF RIGHT SINGULAR VECTORS. */
/*                   V IS NOT REFERENCED IF JOB.EQ.0.  IF P.LE.N, */
/*                   THEN V MAY BE IDENTIFIED WITH X IN THE */
/*                   SUBROUTINE CALL. */

/*         INFO      INTEGER. */
/*                   THE SINGULAR VALUES (AND THEIR CORRESPONDING */
/*                   SINGULAR VECTORS) S(INFO+1),S(INFO+2),...,S(M) */
/*                   ARE CORRECT (HERE M=MIN(N,P)).  THUS IF */
/*                   INFO.EQ.0, ALL THE SINGULAR VALUES AND THEIR */
/*                   VECTORS ARE CORRECT.  IN ANY EVENT, THE MATRIX */
/*                   B = TRANS(U)*X*V IS THE BIDIAGONAL MATRIX */
/*                   WITH THE ELEMENTS OF S ON ITS DIAGONAL AND THE */
/*                   ELEMENTS OF E ON ITS SUPER-DIAGONAL (TRANS(U) */
/*                   IS THE TRANSPOSE OF U).  THUS THE SINGULAR */
/*                   VALUES OF X AND B ARE THE SAME. */

/*     LINPACK. THIS VERSION DATED 03/19/79 . */
/*     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB. */

/*     DSVDC USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS. */

/*     EXTERNAL DROT */
/*     BLAS DAXPY,DDOT,DSCAL,DSWAP,DNRM2,DROTG */
/*     FORTRAN DABS,DMAX1,MAX0,MIN0,MOD,DSQRT */

/*     INTERNAL VARIABLES */



/*     SET THE MAXIMUM NUMBER OF ITERATIONS. */

    /* Parameter adjustments */
    x_dim1 = *ldx;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    --s;
    --e;
    u_dim1 = *ldu;
    u_offset = u_dim1 + 1;
    u -= u_offset;
    v_dim1 = *ldv;
    v_offset = v_dim1 + 1;
    v -= v_offset;
    --work;

    /* Function Body */
    maxit = 30;

/*     DETERMINE WHAT IS TO BE COMPUTED. */

    wantu = FALSE_;
    wantv = FALSE_;
    jobu = *job % 100 / 10;
    ncu = *n;
    if (jobu > 1) {
	ncu = min(*n,*p);
    }
    if (jobu != 0) {
	wantu = TRUE_;
    }
    if (*job % 10 != 0) {
	wantv = TRUE_;
    }

/*     REDUCE X TO BIDIAGONAL FORM, STORING THE DIAGONAL ELEMENTS */
/*     IN S AND THE SUPER-DIAGONAL ELEMENTS IN E. */

    *info = 0;
/* Computing MIN */
    i__1 = *n - 1;
    nct = min(i__1,*p);
/* Computing MAX */
/* Computing MIN */
    i__3 = *p - 2;
    i__1 = 0, i__2 = min(i__3,*n);
    nrt = max(i__1,i__2);
    lu = max(nct,nrt);
    if (lu < 1) {
	goto L170;
    }
    i__1 = lu;
    for (l = 1; l <= i__1; ++l) {
	lp1 = l + 1;
	if (l > nct) {
	    goto L20;
	}

/*           COMPUTE THE TRANSFORMATION FOR THE L-TH COLUMN AND */
/*           PLACE THE L-TH DIAGONAL IN S(L). */

	i__2 = *n - l + 1;
	s[l] = dnrm2_(&i__2, &x[l + l * x_dim1], &c__1);
	if (s[l] == 0.) {
	    goto L10;
	}
	if (x[l + l * x_dim1] != 0.) {
	    s[l] = d_sign(&s[l], &x[l + l * x_dim1]);
	}
	i__2 = *n - l + 1;
	d__1 = 1. / s[l];
	dscal_(&i__2, &d__1, &x[l + l * x_dim1], &c__1);
	x[l + l * x_dim1] += 1.;
L10:
	s[l] = -s[l];
L20:
	if (*p < lp1) {
	    goto L50;
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    if (l > nct) {
		goto L30;
	    }
	    if (s[l] == 0.) {
		goto L30;
	    }

/*              APPLY THE TRANSFORMATION. */

	    i__3 = *n - l + 1;
	    t = -ddot_(&i__3, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &
		    c__1) / x[l + l * x_dim1];
	    i__3 = *n - l + 1;
	    daxpy_(&i__3, &t, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &
		    c__1);
L30:

/*           PLACE THE L-TH ROW OF X INTO  E FOR THE */
/*           SUBSEQUENT CALCULATION OF THE ROW TRANSFORMATION. */

	    e[j] = x[l + j * x_dim1];
/* L40: */
	}
L50:
	if (! wantu || l > nct) {
	    goto L70;
	}

/*           PLACE THE TRANSFORMATION IN U FOR SUBSEQUENT BACK */
/*           MULTIPLICATION. */

	i__2 = *n;
	for (i = l; i <= i__2; ++i) {
	    u[i + l * u_dim1] = x[i + l * x_dim1];
/* L60: */
	}
L70:
	if (l > nrt) {
	    goto L150;
	}

/*           COMPUTE THE L-TH ROW TRANSFORMATION AND PLACE THE */
/*           L-TH SUPER-DIAGONAL IN E(L). */

	i__2 = *p - l;
	e[l] = dnrm2_(&i__2, &e[lp1], &c__1);
	if (e[l] == 0.) {
	    goto L80;
	}
	if (e[lp1] != 0.) {
	    e[l] = d_sign(&e[l], &e[lp1]);
	}
	i__2 = *p - l;
	d__1 = 1. / e[l];
	dscal_(&i__2, &d__1, &e[lp1], &c__1);
	e[lp1] += 1.;
L80:
	e[l] = -e[l];
	if (lp1 > *n || e[l] == 0.) {
	    goto L120;
	}

/*              APPLY THE TRANSFORMATION. */

	i__2 = *n;
	for (i = lp1; i <= i__2; ++i) {
	    work[i] = 0.;
/* L90: */
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *n - l;
	    daxpy_(&i__3, &e[j], &x[lp1 + j * x_dim1], &c__1, &work[lp1], &
		    c__1);
/* L100: */
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *n - l;
	    d__1 = -e[j] / e[lp1];
	    daxpy_(&i__3, &d__1, &work[lp1], &c__1, &x[lp1 + j * x_dim1], &
		    c__1);
/* L110: */
	}
L120:
	if (! wantv) {
	    goto L140;
	}

/*              PLACE THE TRANSFORMATION IN V FOR SUBSEQUENT */
/*              BACK MULTIPLICATION. */

	i__2 = *p;
	for (i = lp1; i <= i__2; ++i) {
	    v[i + l * v_dim1] = e[i];
/* L130: */
	}
L140:
L150:
/* L160: */
	;
    }
L170:

/*     SET UP THE FINAL BIDIAGONAL MATRIX OR ORDER M. */

/* Computing MIN */
    i__1 = *p, i__2 = *n + 1;
    m = min(i__1,i__2);
    nctp1 = nct + 1;
    nrtp1 = nrt + 1;
    if (nct < *p) {
	s[nctp1] = x[nctp1 + nctp1 * x_dim1];
    }
    if (*n < m) {
	s[m] = 0.;
    }
    if (nrtp1 < m) {
	e[nrtp1] = x[nrtp1 + m * x_dim1];
    }
    e[m] = 0.;

/*     IF REQUIRED, GENERATE U. */

    if (! wantu) {
	goto L300;
    }
    if (ncu < nctp1) {
	goto L200;
    }
    i__1 = ncu;
    for (j = nctp1; j <= i__1; ++j) {
	i__2 = *n;
	for (i = 1; i <= i__2; ++i) {
	    u[i + j * u_dim1] = 0.;
/* L180: */
	}
	u[j + j * u_dim1] = 1.;
/* L190: */
    }
L200:
    if (nct < 1) {
	goto L290;
    }
    i__1 = nct;
    for (ll = 1; ll <= i__1; ++ll) {
	l = nct - ll + 1;
	if (s[l] == 0.) {
	    goto L250;
	}
	lp1 = l + 1;
	if (ncu < lp1) {
	    goto L220;
	}
	i__2 = ncu;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *n - l + 1;
	    t = -ddot_(&i__3, &u[l + l * u_dim1], &c__1, &u[l + j * u_dim1], &
		    c__1) / u[l + l * u_dim1];
	    i__3 = *n - l + 1;
	    daxpy_(&i__3, &t, &u[l + l * u_dim1], &c__1, &u[l + j * u_dim1], &
		    c__1);
/* L210: */
	}
L220:
	i__2 = *n - l + 1;
	dscal_(&i__2, &c_b90, &u[l + l * u_dim1], &c__1);
	u[l + l * u_dim1] += 1.;
	lm1 = l - 1;
	if (lm1 < 1) {
	    goto L240;
	}
	i__2 = lm1;
	for (i = 1; i <= i__2; ++i) {
	    u[i + l * u_dim1] = 0.;
/* L230: */
	}
L240:
	goto L270;
L250:
	i__2 = *n;
	for (i = 1; i <= i__2; ++i) {
	    u[i + l * u_dim1] = 0.;
/* L260: */
	}
	u[l + l * u_dim1] = 1.;
L270:
/* L280: */
	;
    }
L290:
L300:

/*     IF IT IS REQUIRED, GENERATE V. */

    if (! wantv) {
	goto L350;
    }
    i__1 = *p;
    for (ll = 1; ll <= i__1; ++ll) {
	l = *p - ll + 1;
	lp1 = l + 1;
	if (l > nrt) {
	    goto L320;
	}
	if (e[l] == 0.) {
	    goto L320;
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *p - l;
	    t = -ddot_(&i__3, &v[lp1 + l * v_dim1], &c__1, &v[lp1 + j * 
		    v_dim1], &c__1) / v[lp1 + l * v_dim1];
	    i__3 = *p - l;
	    daxpy_(&i__3, &t, &v[lp1 + l * v_dim1], &c__1, &v[lp1 + j * 
		    v_dim1], &c__1);
/* L310: */
	}
L320:
	i__2 = *p;
	for (i = 1; i <= i__2; ++i) {
	    v[i + l * v_dim1] = 0.;
/* L330: */
	}
	v[l + l * v_dim1] = 1.;
/* L340: */
    }
L350:

/*     MAIN ITERATION LOOP FOR THE SINGULAR VALUES. */

    mm = m;
    iter = 0;
L360:

/*        QUIT IF ALL THE SINGULAR VALUES HAVE BEEN FOUND. */

/*     ...EXIT */
    if (m == 0) {
	goto L620;
    }

/*        IF TOO MANY ITERATIONS HAVE BEEN PERFORMED, SET */
/*        FLAG AND RETURN. */

    if (iter < maxit) {
	goto L370;
    }
    *info = m;
/*     ......EXIT */
    goto L620;
L370:

/*        THIS SECTION OF THE PROGRAM INSPECTS FOR */
/*        NEGLIGIBLE ELEMENTS IN THE S AND E ARRAYS.  ON */
/*        COMPLETION THE VARIABLES KASE AND L ARE SET AS FOLLOWS. */

/*           KASE = 1     IF S(M) AND E(L-1) ARE NEGLIGIBLE AND L.LT.M */
/*           KASE = 2     IF S(L) IS NEGLIGIBLE AND L.LT.M */
/*           KASE = 3     IF E(L-1) IS NEGLIGIBLE, L.LT.M, AND */
/*                        S(L), ..., S(M) ARE NOT NEGLIGIBLE (QR STEP). */
/*           KASE = 4     IF E(M-1) IS NEGLIGIBLE (CONVERGENCE). */

    i__1 = m;
    for (ll = 1; ll <= i__1; ++ll) {
	l = m - ll;
/*        ...EXIT */
	if (l == 0) {
	    goto L400;
	}
	test = (d__1 = s[l], fabs(d__1)) + (d__2 = s[l + 1], fabs(d__2));
	ztest = test + (d__1 = e[l], fabs(d__1));
	if (ztest != test) {
	    goto L380;
	}
	e[l] = 0.;
/*        ......EXIT */
	goto L400;
L380:
/* L390: */
	;
    }
L400:
    if (l != m - 1) {
	goto L410;
    }
    kase = 4;
    goto L480;
L410:
    lp1 = l + 1;
    mp1 = m + 1;
    i__1 = mp1;
    for (lls = lp1; lls <= i__1; ++lls) {
	ls = m - lls + lp1;
/*           ...EXIT */
	if (ls == l) {
	    goto L440;
	}
	test = 0.;
	if (ls != m) {
	    test += (d__1 = e[ls], fabs(d__1));
	}
	if (ls != l + 1) {
	    test += (d__1 = e[ls - 1], fabs(d__1));
	}
	ztest = test + (d__1 = s[ls], fabs(d__1));
	if (ztest != test) {
	    goto L420;
	}
	s[ls] = 0.;
/*           ......EXIT */
	goto L440;
L420:
/* L430: */
	;
    }
L440:
    if (ls != l) {
	goto L450;
    }
    kase = 3;
    goto L470;
L450:
    if (ls != m) {
	goto L460;
    }
    kase = 1;
    goto L470;
L460:
    kase = 2;
    l = ls;
L470:
L480:
    ++l;

/*        PERFORM THE TASK INDICATED BY KASE. */

    switch ((int)kase) {
	case 1:  goto L490;
	case 2:  goto L520;
	case 3:  goto L540;
	case 4:  goto L570;
    }

/*        DEFLATE NEGLIGIBLE S(M). */

L490:
    mm1 = m - 1;
    f = e[m - 1];
    e[m - 1] = 0.;
    i__1 = mm1;
    for (kk = l; kk <= i__1; ++kk) {
	k = mm1 - kk + l;
	t1 = s[k];
	drotg_(&t1, &f, &cs, &sn);
	s[k] = t1;
	if (k == l) {
	    goto L500;
	}
	f = -sn * e[k - 1];
	e[k - 1] = cs * e[k - 1];
L500:
	if (wantv) {
	    drot_(p, &v[k * v_dim1 + 1], &c__1, &v[m * v_dim1 + 1], &c__1, &
		    cs, &sn);
	}
/* L510: */
    }
    goto L610;

/*        SPLIT AT NEGLIGIBLE S(L). */

L520:
    f = e[l - 1];
    e[l - 1] = 0.;
    i__1 = m;
    for (k = l; k <= i__1; ++k) {
	t1 = s[k];
	drotg_(&t1, &f, &cs, &sn);
	s[k] = t1;
	f = -sn * e[k];
	e[k] = cs * e[k];
	if (wantu) {
	    drot_(n, &u[k * u_dim1 + 1], &c__1, &u[(l - 1) * u_dim1 + 1], &
		    c__1, &cs, &sn);
	}
/* L530: */
    }
    goto L610;

/*        PERFORM ONE QR STEP. */

L540:

/*           CALCULATE THE SHIFT. */

/* Computing MAX */
    d__6 = (d__1 = s[m], fabs(d__1)), d__7 = (d__2 = s[m - 1], fabs(d__2)), 
	    d__6 = max(d__6,d__7), d__7 = (d__3 = e[m - 1], fabs(d__3)), d__6 =
	     max(d__6,d__7), d__7 = (d__4 = s[l], fabs(d__4)), d__6 = max(d__6,
	    d__7), d__7 = (d__5 = e[l], fabs(d__5));
    scale = max(d__6,d__7);
    sm = s[m] / scale;
    smm1 = s[m - 1] / scale;
    emm1 = e[m - 1] / scale;
    sl = s[l] / scale;
    el = e[l] / scale;
/* Computing 2nd power */
    d__1 = emm1;
    b = ((smm1 + sm) * (smm1 - sm) + d__1 * d__1) / 2.;
/* Computing 2nd power */
    d__1 = sm * emm1;
    c = d__1 * d__1;
    shift = 0.;
    if (b == 0. && c == 0.) {
	goto L550;
    }
/* Computing 2nd power */
    d__1 = b;
    shift = sqrt(d__1 * d__1 + c);
    if (b < 0.) {
	shift = -shift;
    }
    shift = c / (b + shift);
L550:
    f = (sl + sm) * (sl - sm) - shift;
    g = sl * el;

/*           CHASE ZEROS. */

    mm1 = m - 1;
    i__1 = mm1;
    for (k = l; k <= i__1; ++k) {
	drotg_(&f, &g, &cs, &sn);
	if (k != l) {
	    e[k - 1] = f;
	}
	f = cs * s[k] + sn * e[k];
	e[k] = cs * e[k] - sn * s[k];
	g = sn * s[k + 1];
	s[k + 1] = cs * s[k + 1];
	if (wantv) {
	    drot_(p, &v[k * v_dim1 + 1], &c__1, &v[(k + 1) * v_dim1 + 1], &
		    c__1, &cs, &sn);
	}
	drotg_(&f, &g, &cs, &sn);
	s[k] = f;
	f = cs * e[k] + sn * s[k + 1];
	s[k + 1] = -sn * e[k] + cs * s[k + 1];
	g = sn * e[k + 1];
	e[k + 1] = cs * e[k + 1];
	if (wantu && k < *n) {
	    drot_(n, &u[k * u_dim1 + 1], &c__1, &u[(k + 1) * u_dim1 + 1], &
		    c__1, &cs, &sn);
	}
/* L560: */
    }
    e[m - 1] = f;
    ++iter;
    goto L610;

/*        CONVERGENCE. */

L570:

/*           MAKE THE SINGULAR VALUE  POSITIVE. */

    if (s[l] >= 0.) {
	goto L580;
    }
    s[l] = -s[l];
    if (wantv) {
	dscal_(p, &c_b90, &v[l * v_dim1 + 1], &c__1);
    }
L580:

/*           ORDER THE SINGULAR VALUE. */

L590:
    if (l == mm) {
	goto L600;
    }
/*           ...EXIT */
    if (s[l] >= s[l + 1]) {
	goto L600;
    }
    t = s[l];
    s[l] = s[l + 1];
    s[l + 1] = t;
    if (wantv && l < *p) {
	dswap_(p, &v[l * v_dim1 + 1], &c__1, &v[(l + 1) * v_dim1 + 1], &c__1);
    }
    if (wantu && l < *n) {
	dswap_(n, &u[l * u_dim1 + 1], &c__1, &u[(l + 1) * u_dim1 + 1], &c__1);
    }
    ++l;
    goto L590;
L600:
    iter = 0;
    --m;
L610:
    goto L360;
L620:
    return 0;
} /* dsvdc_ */

/* Subroutine */ int drotg_(da, db, c, s)
doublereal *da, *db, *c, *s;
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal r, scale, z, roe;


/*     CONSTRUCT GIVENS PLANE ROTATION. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    roe = *db;
    if (fabs(*da) > fabs(*db)) {
	roe = *da;
    }
    scale = fabs(*da) + fabs(*db);
    if (scale != 0.) {
	goto L10;
    }
    *c = 1.;
    *s = 0.;
    r = 0.;
    goto L20;
L10:
/* Computing 2nd power */
    d__1 = *da / scale;
/* Computing 2nd power */
    d__2 = *db / scale;
    r = scale * sqrt(d__1 * d__1 + d__2 * d__2);
    r = d_sign(&c_b149, &roe) * r;
    *c = *da / r;
    *s = *db / r;
L20:
    z = 1.;
    if (fabs(*da) > fabs(*db)) {
	z = *s;
    }
    if (fabs(*db) >= fabs(*da) && *c != 0.) {
	z = 1. / *c;
    }
    *da = r;
    *db = z;
    return 0;
} /* drotg_ */

/* Subroutine */ int drot_(n, dx, incx, dy, incy, c, s)
integer *n;
doublereal *dx;
integer *incx;
doublereal *dy;
integer *incy;
doublereal *c, *s;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i;
    static doublereal dtemp;
    static integer ix, iy;


/*     APPLIES A PLANE ROTATION. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL */
/*         TO 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	dtemp = *c * dx[ix] + *s * dy[iy];
	dy[iy] = *c * dy[iy] - *s * dx[ix];
	dx[ix] = dtemp;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*       CODE FOR BOTH INCREMENTS EQUAL TO 1 */

L20:
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	dtemp = *c * dx[i] + *s * dy[i];
	dy[i] = *c * dy[i] - *s * dx[i];
	dx[i] = dtemp;
/* L30: */
    }
    return 0;
} /* drot_ */

