/*
 * int_fjt.c
 *
 * Copyright (C) 1996 Limit Point Systems, Inc.
 *
 * Author: Curtis Janssen <cljanss@ca.sandia.gov>
 * Maintainer: LPS
 *
 * This file is part of the SC Toolkit.
 *
 * The SC Toolkit is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Library General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * The SC Toolkit is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public License
 * along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
 * the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * The U.S. Government is granted a limited license as per AL 91-7.
 */

/* These routines are based on the gamfun program of
 * Trygve Ulf Helgaker (fall 1984)
 * and calculates the incomplete gamma function as
 * described by McMurchie & Davidson, J. Comp. Phys. 26 (1978) 218.
 * The original routine computed the function for maximum j = 20.
 */

/* $Log$
 * Revision 1.6  1996/10/25 23:28:07  etseidl
 * add copyleft notice
 *
 * Revision 1.5  1996/03/23 02:37:46  cljanss
 * Everything can now be configured with autoconf.
 *
 * Revision 1.4  1995/03/17 01:49:32  cljanss
 * Removed -I. and -I$(SRCDIR) from the default include path in
 * GlobalMakefile to avoid name conflicts with system include files.
 * Modified files under src.lib to include all files relative to src.lib.
 * Makefiles under src.bin need to add the -I. and -I$(SRCDIR) back onto
 * INCLUDE and CXXINCLUDE or make other arrangements.
 *
 * Revision 1.3  1994/08/26  22:45:34  etseidl
 * fix a bunch of warnings, get rid of rcs id's, get rid of bread/bwrite and
 * fread/fwrite modules
 *
 * Revision 1.2  1993/12/30  13:32:51  etseidl
 * mostly rcs id stuff
 *
 * Revision 1.4  1992/06/17  22:04:46  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.3  1992/04/16  17:05:31  jannsen
 * Some memory was not free'd.  (from Eric Remy of Stanford)
 *
 * Revision 1.2  1992/03/31  01:22:34  jannsen
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.2  1992/03/18  23:14:35  cljanss
 * convert to NCUBE os v 3.0
 *
 * Revision 1.1  1991/12/14  00:19:36  cljanss
 * Initial revision
 * */
/*
 * Revision 1.4  91/11/20  18:10:11  cljanss
 * Modified to overcome strangeness on the NCUBE.
 * Put in a slight performance improvement.
 * 
 * Revision 1.3  91/10/31  14:46:14  cljanss
 * Removed a division
 * 
 * Revision 1.2  91/09/28  19:26:50  cljanss
 * Switch to new file naming scheme
 * 
 * Revision 1.1  91/06/16  16:40:07  janssen
 * Initial revision
 *  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>

#define ALLOC_FJTTABLE
#include <chemistry/qc/intv2/fjttable.h>

#include <chemistry/qc/intv2/int_fjt.gbl>
#include <chemistry/qc/intv2/int_fjt.lcl>

static int maxj;
static double_matrix_t gtable;
static double *denomarray;
static double wval_infinity;
static int itable_infinity;

 /* Tablesize should always be at least 121. */
#define TABLESIZE 121

/* Tabulate the incomplete gamma function and put in gtable. */
/*
 *     For J = JMAX a power series expansion is used, see for
 *     example Eq.(39) given by V. Saunders in "Computational
 *     Techniques in Quantum Chemistry and Molecular Physics",
 *     Reidel 1975.  For J < JMAX the values are calculated
 *     using downward recursion in J.
 */
GLOBAL_FUNCTION void
int_initialize_fjt(max)
int max;
{
  int i,j;
  double denom,d2jmax1,r2jmax1,wval,d2wval,sum,term,rexpw;

  maxj = max;

  /* Allocate storage for gtable and int_fjttable. */
  if (  allocbn_double_matrix(&gtable,"n1 n2",maxj + 7,TABLESIZE)
      ||allocbn_double_vector(&int_fjttable,"n",maxj+1)) {
    fprintf(stderr,"allocation of gamma table of size");
    fprintf(stderr," %dx%d and fjt table of size %d failed\n",
            maxj + 7, TABLESIZE, maxj);
    exit(1);
    }

  /* Tabulate the gamma function for t(=wval)=0.0. */
  denom = 1.0;
  for (i=0; i<gtable.n1; i++) {
    gtable.d[i][0] = 1.0/denom;
    denom += 2.0;
    }

  /* Tabulate the gamma function from t(=wval)=0.1, to 12.0. */
  d2jmax1 = 2.0*(gtable.n1-1) + 1.0;
  r2jmax1 = 1.0/d2jmax1;
  for (i=1; i<TABLESIZE; i++) {
    wval = 0.1 * i;
    d2wval = 2.0 * wval;
    term = r2jmax1;
    sum = term;
    denom = d2jmax1;
    for (j=2; j<=200; j++) {
      denom = denom + 2.0;
      term = term * d2wval / denom;
      sum = sum + term;
      if (term <= 1.0e-15) break;
      }
    rexpw = exp(-wval);

    /* Fill in the values for the highest j gtable entries (top row). */
    gtable.d[gtable.n1-1][i] = rexpw * sum;

    /* Work down the table filling in the rest of the column. */
    denom = d2jmax1;
    for (j=gtable.n1 - 2; j>=0; j--) {
      denom = denom - 2.0;
      gtable.d[j][i] = (gtable.d[j+1][i]*d2wval + rexpw)/denom;
      }
    }

  /* Form some denominators, so divisions can be eliminated below. */
  denomarray = (double *) malloc(sizeof(double)*(max+1));
  denomarray[0] = 0.0;
  for (i=1; i<=max; i++) {
    denomarray[i] = 1.0/(2*i - 1);
    }

  wval_infinity = 2*max + 37.0;
  itable_infinity = (int) (10 * wval_infinity);

  }

/* This is called when the fjt routines are no longer needed, or
 * before they are reinitialized with a new maxj. */
GLOBAL_FUNCTION void
int_done_fjt()
{
  free_double_matrix(&gtable);
  free_double_vector(&int_fjttable);
  free(denomarray);
  }

/* Using the tabulated incomplete gamma function in gtable, compute
 * the incomplete gamma function for a particular wval for all 0<=j<=J.
 * The result is placed in the global intermediate int_fjttable.
 */
GLOBAL_FUNCTION void
int_fjt(J,wval)
int J;
double wval;
{
  const double sqrpih =  0.886226925452758;
  const double coef2 =  0.5000000000000000;
  const double coef3 = -0.1666666666666667;
  const double coef4 =  0.0416666666666667;
  const double coef5 = -0.0083333333333333;
  const double coef6 =  0.0013888888888889;
  const double gfac30 =  0.4999489092;
  const double gfac31 = -0.2473631686;
  const double gfac32 =  0.321180909;
  const double gfac33 = -0.3811559346;
  const double gfac20 =  0.4998436875;
  const double gfac21 = -0.24249438;
  const double gfac22 =  0.24642845;
  const double gfac10 =  0.499093162;
  const double gfac11 = -0.2152832;
  const double gfac00 = -0.490;

  double wdif, d2wal, rexpw, /* denom, */ gval, factor, rwval, term;
  int i, itable, irange;

  if (J>maxj) {
    fprintf(stderr,"the int_fjt routine has been incorrectly used\n");
    fprintf(stderr,"J = %d but maxj = %d\n",J,maxj);
    exit(1);
    }

  /* Compute an index into the table. */
  /* The test is needed to avoid floating point exceptions for
   * large values of wval. */
  if (wval > wval_infinity) {
    itable = itable_infinity;
    }
  else {
    itable = (int) (10.0 * wval);
    }

  /* If itable is small enough use the table to compute int_fjttable. */
  if (itable < TABLESIZE) {

    wdif = wval - 0.1 * itable;

    /* Compute fjt for J. */
    int_fjttable.d[J] = (((((coef6 * gtable.d[J+6][itable]*wdif
                            + coef5 * gtable.d[J+5][itable])*wdif
                             + coef4 * gtable.d[J+4][itable])*wdif
                              + coef3 * gtable.d[J+3][itable])*wdif
                               + coef2 * gtable.d[J+2][itable])*wdif
                                -  gtable.d[J+1][itable])*wdif
                         +  gtable.d[J][itable];

    /* Compute the rest of the fjt. */
    d2wal = 2.0 * wval;
    rexpw = exp(-wval);
    /* denom = 2*J + 1; */
    for (i=J; i>0; i--) {
      /* denom = denom - 2.0; */
      int_fjttable.d[i-1] = (d2wal*int_fjttable.d[i] + rexpw)*denomarray[i];
      }
    }
  /* If wval <= 2*J + 36.0, use the following formula. */
  else if (itable <= 20*J + 360) {
    rwval = 1.0/wval;
    rexpw = exp(-wval);

    /* Subdivide wval into 6 ranges. */
    irange = itable/30 - 3;
    if (irange == 1) {
      gval = gfac30 + rwval*(gfac31 + rwval*(gfac32 + rwval*gfac33));
      int_fjttable.d[0] = sqrpih*sqrt(rwval) - rexpw*gval*rwval;
      }
    else if (irange == 2) {
      gval = gfac20 + rwval*(gfac21 + rwval*gfac22);
      int_fjttable.d[0] = sqrpih*sqrt(rwval) - rexpw*gval*rwval;
      }
    else if (irange == 3 || irange == 4) {
      gval = gfac10 + rwval*gfac11;
      int_fjttable.d[0] = sqrpih*sqrt(rwval) - rexpw*gval*rwval;
      }
    else if (irange == 5 || irange == 6) {
      gval = gfac00;
      int_fjttable.d[0] = sqrpih*sqrt(rwval) - rexpw*gval*rwval;
      }
    else {
      int_fjttable.d[0] = sqrpih*sqrt(rwval);
      }

    /* Compute the rest of the int_fjttable from int_fjttable.d[0]. */
    factor = 0.5 * rwval;
    term = factor * rexpw;
    for (i=1; i<=J; i++) {
      int_fjttable.d[i] = factor * int_fjttable.d[i-1] - term;
      factor = rwval + factor;
      }
    }
  /* For large values of wval use this algorithm: */
  else {
    rwval = 1.0/wval;
    int_fjttable.d[0] = sqrpih*sqrt(rwval);
    factor = 0.5 * rwval;
    for (i=1; i<=J; i++) {
      int_fjttable.d[i] = factor * int_fjttable.d[i-1];
      factor = rwval + factor;
      }
    }
  /* printf(" %2d %12.8f %4d %12.8f\n",J,wval,itable,int_fjttable.d[0]); */
  }

