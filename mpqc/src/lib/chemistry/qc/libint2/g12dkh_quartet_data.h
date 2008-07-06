//
// g12dkh_quartet_data.h
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

#ifndef _chemistry_qc_libint2_g12dkhquartetdata_h
#define _chemistry_qc_libint2_g12dkhquartetdata_h

#include <util/misc/math.h>
#include <chemistry/qc/libint2/static.h>
#include <chemistry/qc/libint2/libint2_utils.h>


/*--------------------------------------------------------------------------------
  This function computes constants used in OSRR for a given quartet of primitives

  gamma is the exponent of the Gaussian geminal in the integral
 --------------------------------------------------------------------------------*/
inline void G12DKHLibint2::g12dkh_quartet_data_(prim_data *Data, double scale, double gamma)
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

  /* AC */
#if LIBINT2_DEFINED(g12,BD_x)
  Data->BD_x[0] = quartet_info_.B[0] - quartet_info_.D[0];
#endif
#if LIBINT2_DEFINED(g12,BD_y)
  Data->BD_y[0] = quartet_info_.B[1] - quartet_info_.D[1];
#endif
#if LIBINT2_DEFINED(g12,BD_z)
  Data->BD_z[0] = quartet_info_.B[2] - quartet_info_.D[2];
#endif
  /* BD */
#if LIBINT2_DEFINED(g12,BD_x)
  Data->BD_x[0] = quartet_info_.B[0] - quartet_info_.D[0];
#endif
#if LIBINT2_DEFINED(g12,BD_y)
  Data->BD_y[0] = quartet_info_.B[1] - quartet_info_.D[1];
#endif
#if LIBINT2_DEFINED(g12,BD_z)
  Data->BD_z[0] = quartet_info_.B[2] - quartet_info_.D[2];
#endif

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
  double T = rho2 * oorhog * PQ2;

  //
  // (00|0|00), (00|2|00), and (00|4|00) need to start recursion for (ab|0|cd), (ab|2|cd), (ab|4|cd)
  //
  double rorg = rho * oorhog;
  double sqrt_rorg = sqrt(rorg);
  Data->LIBINT_T_SS_K0G12_SS_0[0] = rorg * sqrt_rorg * exp(-gamma*rorg*PQ2) * pfac_normovlp;
  const double ssss_o_rhog = Data->LIBINT_T_SS_K0G12_SS_0[0] * oorhog;
  Data->LIBINT_T_SS_K2G12_SS_0[0] = (1.5 + T) * ssss_o_rhog;
  Data->LIBINT_T_SS_K4G12_SS_0[0] = (15.0/4.0 + 5.0*T + T*T) * ssss_o_rhog * oorhog;

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
