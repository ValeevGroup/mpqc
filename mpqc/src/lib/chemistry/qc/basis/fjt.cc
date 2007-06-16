//
// fjt.cc
//
// Copyright (C) 2001 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
// Maintainer: EV
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
#include <util/misc/math.h>

#include <iostream>
#include <chemistry/qc/libint2/fjt.h>

using namespace std;

/*------------------------------------------------------
  Initialize Taylor_Fm_Eval object (computes incomplete
  gamma function via Taylor interpolation)
 ------------------------------------------------------*/
Taylor_Fjt_Eval::Taylor_Fjt_Eval(unsigned int mmax, double accuracy)
{
  int i, m;
  int T_idx;
  double T, T_new;
  double egamma, func, dfuncdT;
  double term, sum, denom, rel_error;

  cutoff = epsilon;
  /*---------------------------------------
    We are doing Taylor interpolation with
    n=TAYLOR_ORDER terms here:
    error <= delT^n/(n+1)!
   ---------------------------------------*/
  order_interp = TAYLOR_ORDER;
  delT = 2.0*pow(cutoff*fac[order_interp+1],
			1.0/order_interp);
  max_m = mmax + order_interp - 1;

  /*------------------------------------------------
     Check if Taylor_Fm_Eval has been initialized with
     the same mmax before:
     2) yes  - re-initialize again
     3) no - initialize
   ------------------------------------------------*/
  if (grid != NULL || T_crit != NULL) {
    free_Taylor_Fm_Eval();
  }
  
  T_crit = new double[max_m + 1];   /*--- m=0 is included! ---*/
  max_T = 0;
  /*--- Figure out T_crit for each m and put into the T_crit ---*/
  for(m=max_m;m>=0;m--) {
    /*------------------------------------------
      Newton-Raphson method to solve
      T^{m-0.5}*exp(-T) = epsilon*Gamma(m+0.5)
      The solution is the max T for which to do
      the interpolation
     ------------------------------------------*/
    T = -log(epsilon);
    egamma = epsilon*sqrt(M_PI)*df[2*m]/pow(2,m);
    T_new = T;
    do {
      T = T_new;
      func = pow(T,m-0.5) * exp(-T) - egamma;
      dfuncdT = ((m-0.5) * pow(T,m-1.5) - pow(T,m-0.5)) * exp(-T);
      if (dfuncdT >= 0.0) {
	T_new *= 2.5;
	continue;
      }
      T_new = T - func/dfuncdT;
      if ( T_new <= 0.0 ) {
	T_new = T / 2.0;
      }
    } while (fabs(func/egamma) >= SOFT_ZERO);
    T_crit[m] = T_new;
    T_idx = floor(T_new/delT);
    if (T_idx > max_T)
      max_T = T_idx;
  }

  /*-------------------------------------------------------
    Tabulate the gamma function from t=delT to T_crit[m]:
    1) include T=0 though the table is empty for T=0 since
       Fm(0) is simple to compute
    2) modified MacLaurin series converges fastest for
       the largest m -> use it to compute Fmmax(T)
       see JPC 94, 5564 (1990).
    3) then either use the series to compute the rest
       of the row or maybe use downward recursion
   -------------------------------------------------------*/
  grid = block_matrix(max_T+1,max_m+1);
  /*--- do the mmax first ---*/
  for(m=0;m<=max_m;m++)
  for(T_idx = max_T;
      T_idx >= 0;
      T_idx--) {
    T = T_idx*delT;
    denom = (m+0.5);
    term = 0.5*exp(-T)/denom;
    sum = term;
    do {
      denom += 1.0;
      term *= T/denom;
      sum += term;
      rel_error = term/sum;
    } while (rel_error >= cutoff);

    grid[T_idx][m] = sum;
  }

}

Taylor_Fjt_Eval::~Taylor_Fm_Eval()
{
  delete[] T_crit;
  T_crit = NULL;
  free_block(grid);
  grid = NULL;
}

/* Using the tabulated incomplete gamma function in gtable, compute
 * the incomplete gamma function for a particular wval for all 0<=j<=J.
 * The result is placed in the global intermediate int_fjttable.
 */
double *
Taylor_Fjt_Eval::compute_Fjt(double T, unsigned int l)
{

    int m;
  unsigned int T_ind;
  double T_crit, two_T, exp_mT, h, F_m, F_mp1;
  double *F_row;

#define STATIC_OO2NP1
#define STATIC_OON
#include "static.h"

  T_crit = Taylor_Fm_Eval.T_crit[l];
  two_T = 2.0*T;

  /*------------------------
    First compute Fl(T) ...
   ------------------------*/
  if (T > T_crit) {
    /*--- Asymptotic formula ---*/
    F[l] = df[2*l]*sqrt(M_PI/2)/pow(two_T,l+0.5);
  }
  else {
    /*--- Taylor interpolation ---*/
    T_ind = floor(0.5+T/Taylor_Fm_Eval.delT);
    h = T_ind*Taylor_Fm_Eval.delT - T;
    F_row = Taylor_Fm_Eval.grid[T_ind] + l;
    F[l] =          F_row[0] +
	         h*(F_row[1] +
	  oon[2]*h*(F_row[2] +
	  oon[3]*h*(F_row[3] +
	  oon[4]*h*(F_row[4] +
	  oon[5]*h*(F_row[5])))));
  }

  /*------------------------------------
    And then do downward recursion in m
   ------------------------------------*/
  if (l > 0) {
    F_mp1 = F[l];
    exp_mT = exp(-T);
    for(m=l-1;m>=0;m--) {
      F_m = (exp_mT + two_T*F_mp1)*oo2np1[m];
      F[m] = F_m;
      F_mp1 = F_m;
    }
  }

  return Fjt_buffer;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
