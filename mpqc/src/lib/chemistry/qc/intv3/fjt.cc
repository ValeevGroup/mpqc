//
// fjt.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifdef __GNUG__
#pragma implementation
#endif

/* These routines are based on the gamfun program of
 * Trygve Ulf Helgaker (fall 1984)
 * and calculates the incomplete gamma function as
 * described by McMurchie & Davidson, J. Comp. Phys. 26 (1978) 218.
 * The original routine computed the function for maximum j = 20.
 */

#include <stdlib.h>
#include <math.h>

#include <iostream>

#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <chemistry/qc/intv3/fjt.h>

using namespace std;
using namespace sc;

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
FJT::FJT(int max)
{
  int i,j;
  double denom,d2jmax1,r2jmax1,wval,d2wval,sum,term,rexpw;

  maxj = max;

  /* Allocate storage for gtable and int_fjttable. */
  int_fjttable = new double[maxj+1];
  gtable = new double*[ngtable()];
  for (i=0; i<ngtable(); i++) {
      gtable[i] = new double[TABLESIZE];
    }

  /* Tabulate the gamma function for t(=wval)=0.0. */
  denom = 1.0;
  for (i=0; i<ngtable(); i++) {
    gtable[i][0] = 1.0/denom;
    denom += 2.0;
    }

  /* Tabulate the gamma function from t(=wval)=0.1, to 12.0. */
  d2jmax1 = 2.0*(ngtable()-1) + 1.0;
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
    gtable[ngtable()-1][i] = rexpw * sum;

    /* Work down the table filling in the rest of the column. */
    denom = d2jmax1;
    for (j=ngtable() - 2; j>=0; j--) {
      denom = denom - 2.0;
      gtable[j][i] = (gtable[j+1][i]*d2wval + rexpw)/denom;
      }
    }

  /* Form some denominators, so divisions can be eliminated below. */
  denomarray = new double[max+1];
  denomarray[0] = 0.0;
  for (i=1; i<=max; i++) {
    denomarray[i] = 1.0/(2*i - 1);
    }

  wval_infinity = 2*max + 37.0;
  itable_infinity = (int) (10 * wval_infinity);

  }

FJT::~FJT()
{
  delete[] int_fjttable;
  for (int i=0; i<maxj+7; i++) {
      delete[] gtable[i];
    }
  delete[] gtable;
  delete[] denomarray;
  }

/* Using the tabulated incomplete gamma function in gtable, compute
 * the incomplete gamma function for a particular wval for all 0<=j<=J.
 * The result is placed in the global intermediate int_fjttable.
 */
double *
FJT::values(int J,double wval)
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
    ExEnv::errn()
      << scprintf("the int_fjt routine has been incorrectly used")
      << endl;
    ExEnv::errn()
      << scprintf("J = %d but maxj = %d",J,maxj)
      << endl;
    abort();
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
    int_fjttable[J] = (((((coef6 * gtable[J+6][itable]*wdif
                            + coef5 * gtable[J+5][itable])*wdif
                             + coef4 * gtable[J+4][itable])*wdif
                              + coef3 * gtable[J+3][itable])*wdif
                               + coef2 * gtable[J+2][itable])*wdif
                                -  gtable[J+1][itable])*wdif
                         +  gtable[J][itable];

    /* Compute the rest of the fjt. */
    d2wal = 2.0 * wval;
    rexpw = exp(-wval);
    /* denom = 2*J + 1; */
    for (i=J; i>0; i--) {
      /* denom = denom - 2.0; */
      int_fjttable[i-1] = (d2wal*int_fjttable[i] + rexpw)*denomarray[i];
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
      int_fjttable[0] = sqrpih*sqrt(rwval) - rexpw*gval*rwval;
      }
    else if (irange == 2) {
      gval = gfac20 + rwval*(gfac21 + rwval*gfac22);
      int_fjttable[0] = sqrpih*sqrt(rwval) - rexpw*gval*rwval;
      }
    else if (irange == 3 || irange == 4) {
      gval = gfac10 + rwval*gfac11;
      int_fjttable[0] = sqrpih*sqrt(rwval) - rexpw*gval*rwval;
      }
    else if (irange == 5 || irange == 6) {
      gval = gfac00;
      int_fjttable[0] = sqrpih*sqrt(rwval) - rexpw*gval*rwval;
      }
    else {
      int_fjttable[0] = sqrpih*sqrt(rwval);
      }

    /* Compute the rest of the int_fjttable from int_fjttable[0]. */
    factor = 0.5 * rwval;
    term = factor * rexpw;
    for (i=1; i<=J; i++) {
      int_fjttable[i] = factor * int_fjttable[i-1] - term;
      factor = rwval + factor;
      }
    }
  /* For large values of wval use this algorithm: */
  else {
    rwval = 1.0/wval;
    int_fjttable[0] = sqrpih*sqrt(rwval);
    factor = 0.5 * rwval;
    for (i=1; i<=J; i++) {
      int_fjttable[i] = factor * int_fjttable[i-1];
      factor = rwval + factor;
      }
    }
  /* printf(" %2d %12.8f %4d %12.8f\n",J,wval,itable,int_fjttable[0]); */

  return int_fjttable;
  }

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
