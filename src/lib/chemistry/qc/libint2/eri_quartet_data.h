//
// eri_quartet_data.h
//
// Copyright (C) 2001 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
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

#include <util/misc/math.h>
#include <chemistry/qc/libint2/static.h>
#include <chemistry/qc/libint2/libint2_utils.h>

#ifndef _chemistry_qc_libint2_eriquartetdata_h
#define _chemistry_qc_libint2_eriquartetdata_h

namespace sc {

/*--------------------------------------------------------------------------------
  This function computes constants used in OSRR for a given quartet of primitives
 --------------------------------------------------------------------------------*/
inline void EriLibint2::eri_quartet_data_(prim_data *Data, double scale)
{
#define STATIC_OO2NP1
#include "static.h"

  /*----------------
    Local variables
   ----------------*/
  double P[3], Q[3], PQ[3], W[3];
  double small_T = 1E-15;       /*--- Use only one term in Taylor expansion of Fj(T) if T < small_T ---*/

  int p1 = quartet_info_.p1;
  int p2 = quartet_info_.p2;
  int p3 = quartet_info_.p3;
  int p4 = quartet_info_.p4;
  
  double a1 = int_shell1_->exponent(quartet_info_.p1);
  double a2 = int_shell2_->exponent(quartet_info_.p2);
  double a3 = int_shell3_->exponent(quartet_info_.p3);
  double a4 = int_shell4_->exponent(quartet_info_.p4);

  prim_pair_t* pair12;
  prim_pair_t* pair34;
  if (!quartet_info_.p13p24) {
    pair12 = quartet_info_.shell_pair12->prim_pair(*quartet_info_.op1,*quartet_info_.op2);
    pair34 = quartet_info_.shell_pair34->prim_pair(*quartet_info_.op3,*quartet_info_.op4);
  }
  else {
    pair12 = quartet_info_.shell_pair34->prim_pair(*quartet_info_.op3,*quartet_info_.op4);
    pair34 = quartet_info_.shell_pair12->prim_pair(*quartet_info_.op1,*quartet_info_.op2);
  }
  
  double zeta = pair12->gamma;
  double eta = pair34->gamma;
  double ooz = 1.0/zeta;
  double ooe = 1.0/eta;
  double ooze = 1.0/(zeta+eta);
  Data->roz[0] = eta*ooze;
  double rho = zeta*Data->roz[0];
  
  double pfac_norm = int_shell1_->coefficient_unnorm(quartet_info_.gc1,p1)*
  int_shell2_->coefficient_unnorm(quartet_info_.gc2,p2)*
  int_shell3_->coefficient_unnorm(quartet_info_.gc3,p3)*
  int_shell4_->coefficient_unnorm(quartet_info_.gc4,p4);
  double pfac = 2.0*sqrt(rho*M_1_PI)*scale*pair12->ovlp*pair34->ovlp*pfac_norm;

  P[0] = pair12->P[0];
  P[1] = pair12->P[1];
  P[2] = pair12->P[2];
  Q[0] = pair34->P[0];
  Q[1] = pair34->P[1];
  Q[2] = pair34->P[2];
  PQ[0] = P[0] - Q[0];
  PQ[1] = P[1] - Q[1];
  PQ[2] = P[2] - Q[2];
  double PQ2 = PQ[0]*PQ[0];
  PQ2 += PQ[1]*PQ[1];
  PQ2 += PQ[2]*PQ[2];
  double T = rho*PQ2;

  if (!quartet_info_.am) {
    double *fjttable = Fm_Eval_->values(0,T);
    Data->LIBINT_T_SS_EREP_SS(0)[0] = fjttable[0]*pfac;
  }
  else {
    Data->oo2ze[0] = 0.5*ooze;
    Data->roe[0] = zeta*ooze;
    Data->oo2z[0] = 0.5/zeta;
    Data->oo2e[0] = 0.5/eta;
    W[0] = (zeta*P[0] + eta*Q[0])*ooze;
    W[1] = (zeta*P[1] + eta*Q[1])*ooze;
    W[2] = (zeta*P[2] + eta*Q[2])*ooze;
    
    if(T < small_T){
      assign_FjT(Data,quartet_info_.am,oo2np1,pfac);
    }
    else {
      double *fjttable = Fm_Eval_->values(quartet_info_.am,T);
      assign_FjT(Data,quartet_info_.am,fjttable,pfac);
    }

    /* PA */
    Data->PA_x[0] = P[0] - quartet_info_.A[0];
    Data->PA_y[0] = P[1] - quartet_info_.A[1];
    Data->PA_z[0] = P[2] - quartet_info_.A[2];
    /* QC */
    Data->QC_x[0] = Q[0] - quartet_info_.C[0];
    Data->QC_y[0] = Q[1] - quartet_info_.C[1];
    Data->QC_z[0] = Q[2] - quartet_info_.C[2];
    /* WP */
    Data->WP_x[0] = W[0] - P[0];
    Data->WP_y[0] = W[1] - P[1];
    Data->WP_z[0] = W[2] - P[2];
    /* WQ */
    Data->WQ_x[0] = W[0] - Q[0];
    Data->WQ_y[0] = W[1] - Q[1];
    Data->WQ_z[0] = W[2] - Q[2];
  }

  return;
}

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
