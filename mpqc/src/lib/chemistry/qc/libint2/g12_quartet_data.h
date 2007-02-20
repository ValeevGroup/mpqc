//
// g12_quartet_data.h
//
// Copyright (C) 2005 Edward Valeev
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

#include <util/misc/math.h>
#include <chemistry/qc/libint2/static.h>
#include <chemistry/qc/libint2/libint2_utils.h>

#ifndef _chemistry_qc_libint2_g12quartetdata_h
#define _chemistry_qc_libint2_g12quartetdata_h

/*--------------------------------------------------------------------------------
  This function computes constants used in OSRR for a given quartet of primitives

  gamma is the exponent of the Gaussian geminal in the integral
 --------------------------------------------------------------------------------*/
inline void G12Libint2::g12_quartet_data_(prim_data *Data, double scale, double gamma, bool eri_only)
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

  //
  // Prefactors for [T_i,g12] integrals
  //
  Data->zeta_A[0] = a1;
  Data->zeta_B[0] = a2;
  Data->zeta_C[0] = a3;
  Data->zeta_D[0] = a4;
  Data->zeta_A_2[0] = a1*a1;
  Data->zeta_B_2[0] = a2*a2;
  Data->zeta_C_2[0] = a3*a3;
  Data->zeta_D_2[0] = a4*a4;

  //
  // prefactors for (ab|-1|cd) are same as for OSRR, only (00|-1|00)^m are different
  //
  double zeta = pair12->gamma;
  double eta = pair34->gamma;
  double ooz = 1.0/zeta;
  double ooe = 1.0/eta;
  double ooze = 1.0/(zeta+eta);
  Data->roz[0] = eta*ooze;
  double rho = zeta*Data->roz[0];
  double rhog = rho + gamma;
  double oorhog = 1.0/rhog;
  double rho2 = rho*rho;

  P[0] = pair12->P[0];
  P[1] = pair12->P[1];
  P[2] = pair12->P[2];
  Q[0] = pair34->P[0];
  Q[1] = pair34->P[1];
  Q[2] = pair34->P[2];

  Data->oo2ze[0] = 0.5*ooze;
  Data->roe[0] = zeta*ooze;
  Data->oo2z[0] = 0.5 * ooz;
  Data->oo2e[0] = 0.5 * ooe;
  W[0] = (zeta*P[0] + eta*Q[0])*ooze;
  W[1] = (zeta*P[1] + eta*Q[1])*ooze;
  W[2] = (zeta*P[2] + eta*Q[2])*ooze;

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

  PQ[0] = P[0] - Q[0];
  PQ[1] = P[1] - Q[1];
  PQ[2] = P[2] - Q[2];
  double PQ2 = PQ[0]*PQ[0];
  PQ2 += PQ[1]*PQ[1];
  PQ2 += PQ[2]*PQ[2];

  const double pfac_norm = int_shell1_->coefficient_unnorm(quartet_info_.gc1,p1)*
                           int_shell2_->coefficient_unnorm(quartet_info_.gc2,p2)*
                           int_shell3_->coefficient_unnorm(quartet_info_.gc3,p3)*
                           int_shell4_->coefficient_unnorm(quartet_info_.gc4,p4);
  const double pfac_normovlp = pfac_norm * pair12->ovlp * pair34->ovlp * scale;

  if (eri_only) {
      double T = rho*PQ2;
      double pfac = 2.0*sqrt(rho*M_1_PI)*pfac_normovlp;
      if(T < small_T){
          assign_FjT(Data,quartet_info_.am,oo2np1,pfac);
        }
      else {
          double *fjttable = Fm_Eval_->values(quartet_info_.am,T);
          assign_FjT(Data,quartet_info_.am,fjttable,pfac);
        }
      return;
    }

  // else, if need other integrals
  double T = rho2 * oorhog * PQ2;

  //
  // (00|0|00) and (00|2|00) need to start recursion for (ab|0|cd), (ab|2|cd), and [T_i,G12]
  //
  double rorg = rho * oorhog;
  double sqrt_rorg = sqrt(rorg);
  Data->LIBINT_T_SS_K0G12_SS_0[0] = rorg * sqrt_rorg * exp(-gamma*rorg*PQ2) * pfac_normovlp;
  Data->LIBINT_T_SS_K2G12_SS_0[0] = (1.5 + T) * Data->LIBINT_T_SS_K0G12_SS_0[0] * oorhog;

  //
  // compute (00|-1|00)^m from Fj(x)
  //
  double pfac = 2.0 * sqrt(rhog*M_1_PI) * Data->LIBINT_T_SS_K0G12_SS_0[0];

  const double *F;
  if(T < small_T){ 
    F = oo2np1;
  }
  else {
    F = Fm_Eval_->values(quartet_info_.am,T);
  }

  double ss_m1_ss[4*LIBINT2_MAX_AM_R12kG12+1];
  double g_i[4*LIBINT2_MAX_AM_R12kG12+1];
  double r_i[4*LIBINT2_MAX_AM_R12kG12+1];
  double oorhog_i[4*LIBINT2_MAX_AM_R12kG12+1];
  g_i[0] = 1.0;
  r_i[0] = 1.0;
  oorhog_i[0] = 1.0;
  for(int i=1; i<=quartet_info_.am; i++) {
      g_i[i] = g_i[i-1] * gamma;
      r_i[i] = r_i[i-1] * rho;
      oorhog_i[i] = oorhog_i[i-1] * oorhog;
    }
  for(int m=0; m<=quartet_info_.am; m++) {
      double ssss = 0.0;
      for(int k=0; k<=m; k++) {
          ssss += ExpMath_.bc[m][k] * r_i[k] * g_i[m-k] * F[k];
        }
      ss_m1_ss[m] = ssss * oorhog_i[m];
    }
  
  assign_ss_r12m1g12_ss(Data,quartet_info_.am,ss_m1_ss,pfac);


  //
  // prefactors for (a0|k|c0) (k!=-1) VRR
  //
  {
    double u0 = 0.5/(zeta*eta + gamma*(zeta+eta));

    {
      double t00 = a2*(eta + gamma);
      double t01 = gamma*a4;
      double t02 = gamma*eta;
      double T[3];
      for(int w=0;w<3; w++) {
          T[w] = -2.0 * u0 * (t00*(quartet_info_.A[w]-quartet_info_.B[w]) +
                              t01*(quartet_info_.C[w]-quartet_info_.D[w]) +
                              t02*(quartet_info_.A[w]-quartet_info_.C[w]));
        }
      Data->R12kG12_pfac0_0_x[0] = T[0];
      Data->R12kG12_pfac0_0_y[0] = T[1];
      Data->R12kG12_pfac0_0_z[0] = T[2];
    }
    {
      double t00 = a4*(zeta + gamma);
      double t01 = gamma*a2;
      double t02 = gamma*zeta;
      double T[3];
      for(int w=0;w<3; w++) {
          T[w] = -2.0 * u0 * (t00*(quartet_info_.C[w]-quartet_info_.D[w]) +
                              t01*(quartet_info_.A[w]-quartet_info_.B[w]) +
                              t02*(quartet_info_.C[w]-quartet_info_.A[w]));
        }
      Data->R12kG12_pfac0_1_x[0] = T[0];
      Data->R12kG12_pfac0_1_y[0] = T[1];
      Data->R12kG12_pfac0_1_z[0] = T[2];
    }
    {
      Data->R12kG12_pfac1_0[0] = u0 * (eta + gamma);
      Data->R12kG12_pfac1_1[0] = u0 * (zeta + gamma);
    }
    {
      Data->R12kG12_pfac2[0] = u0 * gamma;
    }
    {
      Data->R12kG12_pfac3_0[0] = eta*u0;
      Data->R12kG12_pfac3_1[0] = zeta*u0;
    }
    {
      double T[3];
      for(int w=0;w<3; w++) {
          T[w] = quartet_info_.A[w]-quartet_info_.C[w];
        }
      Data->R12kG12_pfac4_0_x[0] = T[0];
      Data->R12kG12_pfac4_0_y[0] = T[1];
      Data->R12kG12_pfac4_0_z[0] = T[2];
      Data->R12kG12_pfac4_1_x[0] = -T[0];
      Data->R12kG12_pfac4_1_y[0] = -T[1];
      Data->R12kG12_pfac4_1_z[0] = -T[2];
    }
  }

  return;
}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
