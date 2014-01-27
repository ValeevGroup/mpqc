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
#if LIBINT2_DEFINED(eri,roz)
  Data->roz[0] = eta*ooze;
  double rho = zeta*Data->roz[0];
#else
  double rho = zeta * eta * ooze;
#endif
  
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
    Fm_Eval_.eval(Fm_table_, T, 0);
    Data->LIBINT_T_SS_EREP_SS(0)[0] = Fm_table_[0]*pfac;
  }
  else {
#if LIBINT2_DEFINED(eri,oo2ze)
    Data->oo2ze[0] = 0.5*ooze;
#endif
#if LIBINT2_DEFINED(eri,roe)
    Data->roe[0] = zeta*ooze;
#endif
#if LIBINT2_DEFINED(eri,oo2z)
    Data->oo2z[0] = 0.5/zeta;
#endif
#if LIBINT2_DEFINED(eri,oo2e)
    Data->oo2e[0] = 0.5/eta;
#endif
    W[0] = (zeta*P[0] + eta*Q[0])*ooze;
    W[1] = (zeta*P[1] + eta*Q[1])*ooze;
    W[2] = (zeta*P[2] + eta*Q[2])*ooze;
    
    if(T < small_T){
      assign_FjT(Data,quartet_info_.am,oo2np1,pfac);
    }
    else {
      Fm_Eval_.eval(Fm_table_, T, quartet_info_.am);
      assign_FjT(Data,quartet_info_.am,Fm_table_,pfac);
    }

    /* PA */
#if LIBINT2_DEFINED(eri,PA_x)
    Data->PA_x[0] = P[0] - quartet_info_.A[0];
#endif
#if LIBINT2_DEFINED(eri,PA_y)
    Data->PA_y[0] = P[1] - quartet_info_.A[1];
#endif
#if LIBINT2_DEFINED(eri,PA_z)
    Data->PA_z[0] = P[2] - quartet_info_.A[2];
#endif
    /* QC */
#if LIBINT2_DEFINED(eri,QC_x)
    Data->QC_x[0] = Q[0] - quartet_info_.C[0];
#endif
#if LIBINT2_DEFINED(eri,QC_y)
    Data->QC_y[0] = Q[1] - quartet_info_.C[1];
#endif
#if LIBINT2_DEFINED(eri,QC_z)
    Data->QC_z[0] = Q[2] - quartet_info_.C[2];
#endif
    /* WP */
#if LIBINT2_DEFINED(eri,WP_x)
    Data->WP_x[0] = W[0] - P[0];
#endif
#if LIBINT2_DEFINED(eri,WP_y)
    Data->WP_y[0] = W[1] - P[1];
#endif
#if LIBINT2_DEFINED(eri,WP_z)
    Data->WP_z[0] = W[2] - P[2];
#endif
    /* WQ */
#if LIBINT2_DEFINED(eri,WQ_x)
    Data->WQ_x[0] = W[0] - Q[0];
#endif
#if LIBINT2_DEFINED(eri,WQ_y)
    Data->WQ_y[0] = W[1] - Q[1];
#endif
#if LIBINT2_DEFINED(eri,WQ_z)
    Data->WQ_z[0] = W[2] - Q[2];
#endif

    // using ITR?
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_0_x)
    Data->TwoPRepITR_pfac0_0_0_x[0] = - (a2*(quartet_info_.A[0]-quartet_info_.B[0]) + a4*(quartet_info_.C[0]-quartet_info_.D[0]))/zeta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_0_y)
    Data->TwoPRepITR_pfac0_0_0_y[0] = - (a2*(quartet_info_.A[1]-quartet_info_.B[1]) + a4*(quartet_info_.C[1]-quartet_info_.D[1]))/zeta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_0_z)
    Data->TwoPRepITR_pfac0_0_0_z[0] = - (a2*(quartet_info_.A[2]-quartet_info_.B[2]) + a4*(quartet_info_.C[2]-quartet_info_.D[2]))/zeta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_0_x)
    Data->TwoPRepITR_pfac0_1_0_x[0] = - (a2*(quartet_info_.A[0]-quartet_info_.B[0]) + a4*(quartet_info_.C[0]-quartet_info_.D[0]))/eta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_0_y)
    Data->TwoPRepITR_pfac0_1_0_y[0] = - (a2*(quartet_info_.A[1]-quartet_info_.B[1]) + a4*(quartet_info_.C[1]-quartet_info_.D[1]))/eta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_0_z)
    Data->TwoPRepITR_pfac0_1_0_z[0] = - (a2*(quartet_info_.A[2]-quartet_info_.B[2]) + a4*(quartet_info_.C[2]-quartet_info_.D[2]))/eta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_1_x)
    Data->TwoPRepITR_pfac0_0_1_x[0] = (a1*(quartet_info_.A[0]-quartet_info_.B[0]) + a3*(quartet_info_.C[0]-quartet_info_.D[0]))/zeta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_1_y)
    Data->TwoPRepITR_pfac0_0_1_y[0] = (a1*(quartet_info_.A[1]-quartet_info_.B[1]) + a3*(quartet_info_.C[1]-quartet_info_.D[1]))/zeta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_1_z)
    Data->TwoPRepITR_pfac0_0_1_z[0] = (a1*(quartet_info_.A[2]-quartet_info_.B[2]) + a3*(quartet_info_.C[2]-quartet_info_.D[2]))/zeta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_1_x)
    Data->TwoPRepITR_pfac0_1_1_x[0] = (a1*(quartet_info_.A[0]-quartet_info_.B[0]) + a3*(quartet_info_.C[0]-quartet_info_.D[0]))/eta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_1_y)
    Data->TwoPRepITR_pfac0_1_1_y[0] = (a1*(quartet_info_.A[1]-quartet_info_.B[1]) + a3*(quartet_info_.C[1]-quartet_info_.D[1]))/eta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_1_z)
    Data->TwoPRepITR_pfac0_1_1_z[0] = (a1*(quartet_info_.A[2]-quartet_info_.B[2]) + a3*(quartet_info_.C[2]-quartet_info_.D[2]))/eta;
#endif
#if LIBINT2_DEFINED(eri,eoz)
    Data->eoz[0] = eta * ooz;
#endif
#if LIBINT2_DEFINED(eri,zoe)
    Data->zoe[0] = zeta * ooe;
#endif
  }

  return;
}

}
#else

#include <cassert>

namespace {
  class Libint2ExpensiveMath {
  public:
    Libint2ExpensiveMath() {
      const unsigned int imax = 4*LIBINT2_MAX_AM_ERI;
      fac[0] = 1.0;
      for(int i=1; i<=imax; i++) {
        fac[i] = i*fac[i-1];
      }
      for(int i=0; i<=imax; i++) {
        for(int j=0; j<=i; j++) {
          bc[i][j] = fac[i]/(fac[j]*fac[i-j]);
        }
      }
    }
    double fac[4*LIBINT2_MAX_AM_ERI+1];
    double bc[4*LIBINT2_MAX_AM_ERI+1][4*LIBINT2_MAX_AM_ERI+1];
  };
  static Libint2ExpensiveMath libint2_expensive_math;
};

namespace sc {

  // types of integrals that can be computed
  enum OperType {coulomb, f12_coulomb, f12, f12_2, f12_T1_f12};

/*--------------------------------------------------------------------------------
  This function computes constants used in OSRR for a given quartet of primitives

  gamma is the exponents of the Gaussian geminal in the integral

  r12_2_g12_pfac is the prefactor by which the r12^2 * g12 integrals are scaled
  it's necessary to produce [g12,[T1,g12]]
 --------------------------------------------------------------------------------*/
inline void EriLibint2::eri_quartet_data_(prim_data *Data, double scale)
{
#if 0
  OperType otype = f12_T1_f12;
  typedef IntParamsG12::PrimitiveGeminal PrimitiveGeminal;
  typedef IntParamsG12::ContractedGeminal ContractedGeminal;
  ContractedGeminal geminal(3);
  geminal[0] = std::make_pair(1.0, 2.0);
  geminal[1] = std::make_pair(5.0, 3.0);
  geminal[2] = std::make_pair(10.0, 4.0);
  const ContractedGeminal* gbra = &geminal;
  const ContractedGeminal* gket = &geminal;
#endif

#define STATIC_OO2NP1
#include "static.h"

  /*----------------
    Local variables
   ----------------*/
  double P[3], Q[3], PQ[3], W[3];
  const double small_T = 1E-15;       /*--- Use only one term in Taylor expansion of Fj(T) if T < small_T ---*/

  const int p1 = quartet_info_.p1;
  const int p2 = quartet_info_.p2;
  const int p3 = quartet_info_.p3;
  const int p4 = quartet_info_.p4;

  const double a1 = int_shell1_->exponent(quartet_info_.p1);
  const double a2 = int_shell2_->exponent(quartet_info_.p2);
  const double a3 = int_shell3_->exponent(quartet_info_.p3);
  const double a4 = int_shell4_->exponent(quartet_info_.p4);

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
  // prefactors for OSRR do not depend on the integral kernel
  //
  const double zeta = pair12->gamma;
  const double eta = pair34->gamma;
  const double ooz = 1.0/zeta;
  const double ooe = 1.0/eta;
  const double ooze = 1.0/(zeta+eta);
#if LIBINT2_DEFINED(eri,roz)
  Data->roz[0] = eta*ooze;
  double rho = zeta*Data->roz[0];
#else
  double rho = zeta * eta * ooze;
#endif

  const double pfac_norm = int_shell1_->coefficient_unnorm(quartet_info_.gc1,p1)*
                           int_shell2_->coefficient_unnorm(quartet_info_.gc2,p2)*
                           int_shell3_->coefficient_unnorm(quartet_info_.gc3,p3)*
                           int_shell4_->coefficient_unnorm(quartet_info_.gc4,p4);
  const double pfac_normovlp = pfac_norm * pair12->ovlp * pair34->ovlp * scale;


  P[0] = pair12->P[0];
  P[1] = pair12->P[1];
  P[2] = pair12->P[2];
  Q[0] = pair34->P[0];
  Q[1] = pair34->P[1];
  Q[2] = pair34->P[2];
  PQ[0] = P[0] - Q[0];
  PQ[1] = P[1] - Q[1];
  PQ[2] = P[2] - Q[2];
  const double PQ2 = PQ[0]*PQ[0] + PQ[1]*PQ[1] + PQ[2]*PQ[2];
  const double T = rho*PQ2;

  // Coulomb integral
  if (otype == coulomb) {
    double pfac = 2.0*sqrt(rho*M_1_PI)*pfac_normovlp;
    if(T < small_T){
      assign_FjT(Data,quartet_info_.am,oo2np1,pfac);
    }
    else {
      Fm_Eval_.values(Fm_table_, T, quartet_info_.am);
      assign_FjT(Data,quartet_info_.am,Fm_table_,pfac);
    }
  }
  else {

    // this stores (ss|Oper|ss)^(m) integrals
    double ss_oper_ss[4*LIBINT2_MAX_AM_ERI+1];
    std::fill(ss_oper_ss, ss_oper_ss + quartet_info_.am + 1, 0.0);

    // f12_coulomb and f12 integrals
    if (otype == f12 || otype == f12_coulomb) {

      MPQC_ASSERT(gbra != 0);
      const size_t ngbra = gbra->size();
      for(size_t ig=0; ig<ngbra; ++ig) {
        const PrimitiveGeminal& gbra_i = gbra->operator[](ig);
        const double gamma = gbra_i.first;
        const double gcoef = gbra_i.second;
        const double rhog = rho + gamma;
        const double oorhog = 1.0/rhog;
        const double rho2 = rho*rho;

        const double gorg = gamma * oorhog;
        const double rorg = rho * oorhog;
        const double sqrt_rorg = sqrt(rorg);
        /// (ss|g12|ss)
        const double SS_K0G12_SS = rorg * sqrt_rorg * exp(-gorg*T) * pfac_normovlp;

        if (otype == f12_coulomb) {
          const double rorgT = rorg * T;
          double pfac = 2.0 * sqrt(rhog*M_1_PI) * SS_K0G12_SS;

          const double *F;
          if(rorgT < small_T){
            F = oo2np1;
          }
          else {
            Fm_Eval_->values(Fm_table_, rorgT, quartet_info_.am);
            F = Fm_table_;
          }

          double g_i[4*LIBINT2_MAX_AM_ERI+1];
          double r_i[4*LIBINT2_MAX_AM_ERI+1];
          double oorhog_i[4*LIBINT2_MAX_AM_ERI+1];
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
              ssss += libint2_expensive_math.bc[m][k] * r_i[k] * g_i[m-k] * F[k];
            }
            ss_oper_ss[m] += gcoef * pfac * ssss * oorhog_i[m];
          }

        }

        if (otype == f12) {

          double ss_oper_ss_m = SS_K0G12_SS * gcoef;
          ss_oper_ss[0] += ss_oper_ss_m;
          for(int m=1; m<=quartet_info_.am; ++m) {
            ss_oper_ss_m *= gorg;
            ss_oper_ss[m] += ss_oper_ss_m;
          }

        }

      } // loop over gaussian geminals
    } // one correlation factor involved

    // f12_2 and f12_T1_f12 integrals
    if (otype == f12_2 || otype == f12_T1_f12) {

      MPQC_ASSERT(gbra != 0);
      MPQC_ASSERT(gket != 0);
      const size_t ngbra = gbra->size();
      const size_t ngket = gket->size();
      for(size_t igbra=0; igbra<ngbra; ++igbra) {
        const PrimitiveGeminal& gbra_i = gbra->operator[](igbra);
        const double gamma_bra = gbra_i.first;
        const double gcoef_bra = gbra_i.second;
        for(size_t igket=0; igket<ngket; ++igket) {
          const PrimitiveGeminal& gket_i = gket->operator[](igket);
          const double gamma_ket = gket_i.first;
          const double gcoef_ket = gket_i.second;

          const double gamma = gamma_bra + gamma_ket;
          const double gcoef = gcoef_bra * gcoef_ket;

          const double rhog = rho + gamma;
          const double oorhog = 1.0/rhog;
          const double rho2 = rho*rho;

          const double gorg = gamma * oorhog;
          const double rorg = rho * oorhog;
          const double sqrt_rorg = sqrt(rorg);
          /// (ss|g12|ss)
          const double SS_K0G12_SS = rorg * sqrt_rorg * exp(-gorg*T) * pfac_normovlp;
          const double rorgT = rorg * T;
          /// (ss|g12*r12^2|ss)
          const double SS_K2G12_SS_0 = (1.5 + rorgT) * (SS_K0G12_SS * oorhog);
          const double SS_K2G12_SS_m1 = rorg * (SS_K0G12_SS * oorhog);

          if (otype == f12_2) {

            double ss_oper_ss_m = SS_K0G12_SS * gcoef;
            ss_oper_ss[0] += ss_oper_ss_m;
            for(int m=1; m<=quartet_info_.am; ++m) {
              ss_oper_ss_m *= gorg;
              ss_oper_ss[m] += ss_oper_ss_m;
            }

          }

          if (otype == f12_T1_f12) {

            double SS_K2G12_SS_gorg_m = SS_K2G12_SS_0 * (gcoef * 4.0 * gamma_bra * gamma_ket);
            double SS_K2G12_SS_gorg_m1 = SS_K2G12_SS_m1 * (gcoef * 4.0 * gamma_bra * gamma_ket);
            ss_oper_ss[0] += SS_K2G12_SS_gorg_m;
            for(int m=1; m<=quartet_info_.am; ++m) {
              SS_K2G12_SS_gorg_m *= gorg;
              ss_oper_ss[m] += SS_K2G12_SS_gorg_m - m * SS_K2G12_SS_gorg_m1;
              SS_K2G12_SS_gorg_m1 *= gorg;
            }

          }

        }

      } // loop over gaussian geminals
    } // two correlation factors involved

    assign_FjT(Data,quartet_info_.am,ss_oper_ss,1.0);

  }

  // these prefactors only necessary if angular momenta != 0
  if (quartet_info_.am != 0) {
#if LIBINT2_DEFINED(eri,oo2ze)
    Data->oo2ze[0] = 0.5 * ooze;
#endif
#if LIBINT2_DEFINED(eri,roe)
    Data->roe[0] = zeta * ooze;
#endif
#if LIBINT2_DEFINED(eri,oo2z)
    Data->oo2z[0] = 0.5 * ooz;
#endif
#if LIBINT2_DEFINED(eri,oo2e)
    Data->oo2e[0] = 0.5 * ooe;
#endif
    W[0] = (zeta * P[0] + eta * Q[0]) * ooze;
    W[1] = (zeta * P[1] + eta * Q[1]) * ooze;
    W[2] = (zeta * P[2] + eta * Q[2]) * ooze;

    /* PA */
#if LIBINT2_DEFINED(eri,PA_x)
    Data->PA_x[0] = P[0] - quartet_info_.A[0];
#endif
#if LIBINT2_DEFINED(eri,PA_y)
    Data->PA_y[0] = P[1] - quartet_info_.A[1];
#endif
#if LIBINT2_DEFINED(eri,PA_z)
    Data->PA_z[0] = P[2] - quartet_info_.A[2];
#endif
    /* QC */
#if LIBINT2_DEFINED(eri,QC_x)
    Data->QC_x[0] = Q[0] - quartet_info_.C[0];
#endif
#if LIBINT2_DEFINED(eri,QC_y)
    Data->QC_y[0] = Q[1] - quartet_info_.C[1];
#endif
#if LIBINT2_DEFINED(eri,QC_z)
    Data->QC_z[0] = Q[2] - quartet_info_.C[2];
#endif
    /* WP */
#if LIBINT2_DEFINED(eri,WP_x)
    Data->WP_x[0] = W[0] - P[0];
#endif
#if LIBINT2_DEFINED(eri,WP_y)
    Data->WP_y[0] = W[1] - P[1];
#endif
#if LIBINT2_DEFINED(eri,WP_z)
    Data->WP_z[0] = W[2] - P[2];
#endif
    /* WQ */
#if LIBINT2_DEFINED(eri,WQ_x)
    Data->WQ_x[0] = W[0] - Q[0];
#endif
#if LIBINT2_DEFINED(eri,WQ_y)
    Data->WQ_y[0] = W[1] - Q[1];
#endif
#if LIBINT2_DEFINED(eri,WQ_z)
    Data->WQ_z[0] = W[2] - Q[2];
#endif

    // using ITR?
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_0_x)
    Data->TwoPRepITR_pfac0_0_0_x[0] = - (a2*(quartet_info_.A[0]-quartet_info_.B[0]) + a4*(quartet_info_.C[0]-quartet_info_.D[0]))/zeta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_0_y)
    Data->TwoPRepITR_pfac0_0_0_y[0] = - (a2*(quartet_info_.A[1]-quartet_info_.B[1]) + a4*(quartet_info_.C[1]-quartet_info_.D[1]))/zeta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_0_z)
    Data->TwoPRepITR_pfac0_0_0_z[0] = - (a2*(quartet_info_.A[2]-quartet_info_.B[2]) + a4*(quartet_info_.C[2]-quartet_info_.D[2]))/zeta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_0_x)
    Data->TwoPRepITR_pfac0_1_0_x[0] = - (a2*(quartet_info_.A[0]-quartet_info_.B[0]) + a4*(quartet_info_.C[0]-quartet_info_.D[0]))/eta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_0_y)
    Data->TwoPRepITR_pfac0_1_0_y[0] = - (a2*(quartet_info_.A[1]-quartet_info_.B[1]) + a4*(quartet_info_.C[1]-quartet_info_.D[1]))/eta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_0_z)
    Data->TwoPRepITR_pfac0_1_0_z[0] = - (a2*(quartet_info_.A[2]-quartet_info_.B[2]) + a4*(quartet_info_.C[2]-quartet_info_.D[2]))/eta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_1_x)
    Data->TwoPRepITR_pfac0_0_1_x[0] = (a1*(quartet_info_.A[0]-quartet_info_.B[0]) + a3*(quartet_info_.C[0]-quartet_info_.D[0]))/zeta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_1_y)
    Data->TwoPRepITR_pfac0_0_1_y[0] = (a1*(quartet_info_.A[1]-quartet_info_.B[1]) + a3*(quartet_info_.C[1]-quartet_info_.D[1]))/zeta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_1_z)
    Data->TwoPRepITR_pfac0_0_1_z[0] = (a1*(quartet_info_.A[2]-quartet_info_.B[2]) + a3*(quartet_info_.C[2]-quartet_info_.D[2]))/zeta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_1_x)
    Data->TwoPRepITR_pfac0_1_1_x[0] = (a1*(quartet_info_.A[0]-quartet_info_.B[0]) + a3*(quartet_info_.C[0]-quartet_info_.D[0]))/eta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_1_y)
    Data->TwoPRepITR_pfac0_1_1_y[0] = (a1*(quartet_info_.A[1]-quartet_info_.B[1]) + a3*(quartet_info_.C[1]-quartet_info_.D[1]))/eta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_1_z)
    Data->TwoPRepITR_pfac0_1_1_z[0] = (a1*(quartet_info_.A[2]-quartet_info_.B[2]) + a3*(quartet_info_.C[2]-quartet_info_.D[2]))/eta;
#endif
#if LIBINT2_DEFINED(eri,eoz)
    Data->eoz[0] = eta * ooz;
#endif
#if LIBINT2_DEFINED(eri,zoe)
    Data->zoe[0] = zeta * ooe;
#endif

  }

  return;
}

}

#endif // header guard

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
