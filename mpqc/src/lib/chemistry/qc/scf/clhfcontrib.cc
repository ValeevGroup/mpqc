//
// clhfcontrib.cc --- compute the CLHF contribution to the Fock matrix
//
// Based on code:
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
// Maintainer: SNL
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <chemistry/qc/scf/clhfcontrib.h>

#include <iomanip>

#undef DEBUG
#define DEBUG 0

namespace sc {

/////////////////////////////////////////////////////////////////
// CLHFContribution

CLHFContribution::CLHFContribution(
    Ref<GaussianBasisSet> &bs1,
    Ref<GaussianBasisSet> &bs2,
    Ref<GaussianBasisSet> &bs3,
    Ref<GaussianBasisSet> &bs4
    ):
  GenericFockContribution(1,1,bs1,bs2,bs3,bs4)
{
}

CLHFContribution::CLHFContribution(const CLHFContribution &c):
  GenericFockContribution(c)
{
}

CLHFContribution::~CLHFContribution()
{
}

Ref<FockContribution>
CLHFContribution::clone()
{
  return new CLHFContribution(*this);
}

#if DEBUG
static inline void
F_contrib(int I,int J,int K,int L, int i, int j, int k, int l,
          double factor,double integral, const char *type,
          double &F,int fI,int fJ,int fij,int fdim,
          const double &P,int pI,int pJ,int pij,int pdim)
{
  double contrib = factor * integral * P;
  F += contrib;
  if (fabs(contrib) > 1e-6) {
      int fi = fij/fdim, fj = fij%fdim;
      int pi = pij/pdim, pj = pij%pdim;
      std::cout << scprintf("%4.1f * ", factor)
                << "I_" << I << J << K << L
                << "(" << i << j << k << l << ")"
                << " * P_" << pI << pJ
                << "(" << pi << pj << ")"
                << " -> F_" << fI << fJ
                << "(" << fi << fj << ")"
                << " " << type << " " << scprintf("%12.9f",contrib)
                << " val = " << scprintf("%12.9f",F)
                << std::endl;
    }
}
#else
#define \
F_contrib(I,J,K,L,i,j,k,l,\
          factor,integral,type,\
          F,fI,fJ,fij,fdim,\
          P,pI,pJ,pij,pdim)\
  F += factor*integral*P
#endif
    
void
CLHFContribution::contrib_e_J(double factor,
                              int I, int J, int K, int L,
                              int nI, int nJ, int nK, int nL,
                              const double * restrictxx buf)
{
  if (I < J && jmat_symmetric(0)) return;

  double *F_IJ = jmat_block(0, I, J);

  if (K>=L || !pmat_symmetric(0)) {
      const double * restrictxx P_KL = pmat_block(0, K, L);
      int nKL = nK*nL;
      for (int i=0, ijkl=0, ij=0; i<nI; i++) {
          for (int j=0; j<nJ; j++, ij++) {
              double F_IJ_ij = 0.0;
              for (int kl=0; kl<nKL; ijkl++, kl++) {
                  double val = buf[ijkl];
                  F_contrib(I,J,K,L,i,j,kl/nL,kl%nL,factor,val,"J ",
                            F_IJ_ij,I,J,ij,nJ,
                            P_KL[kl],K,L,kl,nL);
                }
              F_IJ[ij] += F_IJ_ij;
            }
        }
    }
  else {
      const double * restrictxx P_LK = pmat_block(0, L, K);
      int nKL = nK*nL;
      for (int i=0, ijkl=0, ij=0; i<nI; i++) {
          for (int j=0; j<nJ; j++, ij++) {
              double F_IJ_ij = 0.0;
              for (int k=0; k<nK; k++) {
                  for (int l=0, lk=k; l<nL; ijkl++, l++, lk+=nK) {
                      double val = buf[ijkl];
                      F_contrib(I,J,K,L,i,j,kl/nL,kl%nL,factor,val,"J ",
                                F_IJ_ij,I,J,ij,nJ,
                                P_LK[lk],L,K,lk,nK);
                    }
                }
              F_IJ[ij] += F_IJ_ij;
            }
        }
    }
}

void
CLHFContribution::contrib_e_K(double factor,
                              int I, int J, int K, int L,
                              int nI, int nJ, int nK, int nL,
                              const double * restrictxx buf)
{
  if (I < K && kmat_symmetric(0)) return;

  double K_factor = factor * -0.5;
  double *F_IK = kmat_block(0, I, K);

  if (J>=L || !pmat_symmetric(0)) {
      const double * restrictxx P_JL = pmat_block(0, J, L);
      for (int i=0, ijkl=0, ik_begin=0; i<nI; i++, ik_begin += nK) {
          for (int j=0, jl_begin=0; j<nJ; j++, jl_begin += nL) {
              for (int k=0, ik=ik_begin; k<nK; k++, ik++) {
                  for (int l=0, jl=jl_begin; l<nL; l++, ijkl++, jl++) {
                      double val = buf[ijkl];
                      F_contrib(I,J,K,L,i,j,k,l,K_factor,val,"K ",
                                F_IK[ik],I,K,ik,nK,
                                P_JL[jl],J,L,jl,nL);
                    }
                }
            }
        }
    }
  else {
      const double * restrictxx P_LJ = pmat_block(0, L, J);
      for (int i=0, ijkl=0, ik_begin=0; i<nI; i++, ik_begin += nK) {
          for (int j=0, lj_begin=0; j<nJ; j++, lj_begin += 1) {
              for (int k=0, ik=ik_begin; k<nK; k++, ik++) {
                  for (int l=0, lj=lj_begin; l<nL; l++, ijkl++, lj+=nJ) {
                      double val = buf[ijkl];
                      F_contrib(I,J,K,L,i,j,k,l,K_factor,val,"K ",
                                F_IK[ik],I,K,ik,nK,
                                P_LJ[lj],L,J,lj,nJ);
                    }
                }
            }
        }
    }
}

void
CLHFContribution::contrib_all_J(double factor,
                                int I, int J, int K, int L,
                                int nI, int nJ, int nK, int nL,
                                const double * restrictxx buf)
{
  double J_factor = factor * 2.0;

  double *F_IJ = jmat_block(0, I, J);
  double *F_KL = jmat_block(0, K, L);
  const double * restrictxx P_IJ = pmat_block(0, I, J);
  const double * restrictxx P_KL = pmat_block(0, K, L);

  int nKL = nK * nL;

  for (int i=0, ijkl=0, ij=0; i<nI; i++) {
      for (int j=0; j<nJ; j++, ij++) {
          double F_IJ_ij = F_IJ[ij];
          double P_IJ_ij = P_IJ[ij];
          for (int kl=0; kl<nKL; kl++, ijkl++) {
              double val = buf[ijkl];
              F_contrib(I,J,K,L,i,j,kl/nL,kl%nL,J_factor,val,"J ",
                        F_IJ_ij,I,J,ij,nJ,
                        P_KL[kl],K,L,kl,nL);
              F_contrib(I,J,K,L,i,j,kl/nL,kl%nL,J_factor,val,"J ",
                        F_KL[kl],K,L,kl,nL,
                        P_IJ_ij,I,J,ij,nJ);
            }
          F_IJ[ij] = F_IJ_ij;
        }
    }
}

void
CLHFContribution::contrib_all_K(double factor,
                                int I, int J, int K, int L,
                                int nI, int nJ, int nK, int nL,
                                const double * restrictxx buf)
{
  double *F_IK = kmat_block(0, I, K);
  double *F_IL = kmat_block(0, I, L);
  const double * restrictxx P_IK = pmat_block(0, I, K);
  const double * restrictxx P_IL = pmat_block(0, I, L);

  // We normally leave out the transposed Fock matrix element and apply the
  // standard contribution of -0.5.  However, if both indices are the same,
  // then we need to multiply by two.
  double IK_factor = -0.5*(I==K?2:1)*factor;
  double JK_factor = -0.5*(J==K?2:1)*factor;
  double IL_factor = -0.5*(I==L?2:1)*factor;
  double JL_factor = -0.5*(J==L?2:1)*factor;

  if (J >= K) {
      double *F_JK = kmat_block(0, J, K);
      double *F_JL = kmat_block(0, J, L);
      const double * restrictxx P_JK = pmat_block(0, J, K);
      const double * restrictxx P_JL = pmat_block(0, J, L);
      for (int i=0,ijkl=0,il_begin=0,ik_begin=0; i<nI;
           i++,il_begin+=nL,ik_begin+=nK) {
          for (int j=0, jk=0, jl_begin=0; j<nJ; j++, jl_begin += nL) {
              for (int k=0,ik=ik_begin; k<nK; k++,ik++,jk++) {
                  for (int l=0, il=il_begin, jl=jl_begin;
                       l<nL;
                       l++, ijkl++, il++, jl++) {
                      double val = buf[ijkl];
                      F_contrib(I,J,K,L,i,j,k,l,IK_factor,val,"K ",
                                F_IK[ik],I,K,ik,nK,
                                P_JL[jl],J,L,jl,nL);
                      F_contrib(I,J,K,L,i,j,k,l,JK_factor,val,"K ",
                                F_JK[jk],J,K,jk,nK,
                                P_IL[il],I,L,il,nL);
                      F_contrib(I,J,K,L,i,j,k,l,IL_factor,val,"K ",
                                F_IL[il],I,L,il,nL,
                                P_JK[jk],J,K,jk,nK);
                      F_contrib(I,J,K,L,i,j,k,l,JL_factor,val,"K ",
                                F_JL[jl],J,L,jl,nL,
                                P_IK[ik],I,K,ik,nK);
                    }
                }
            }
        }
    }
  else if (J >= L) {
      double *F_KJ = kmat_block(0, K, J);
      double *F_JL = kmat_block(0, J, L);
      const double * restrictxx P_KJ = pmat_block(0, K, J);
      const double * restrictxx P_JL = pmat_block(0, J, L);
      for (int i=0,ijkl=0,il_begin=0,ik_begin=0; i<nI;
           i++,il_begin+=nL,ik_begin+=nK) {
          for (int j=0, jl_begin=0; j<nJ; j++, jl_begin += nL) {
              int kj = j;
              for (int k=0,kj=j,ik=ik_begin; k<nK; k++,ik++,kj+=nJ) {
                  for (int l=0, il=il_begin, jl=jl_begin;
                       l<nL;
                       l++, ijkl++, il++, jl++) {
                      double val = buf[ijkl];
                      F_contrib(I,J,K,L,i,j,k,l,IK_factor,val,"K ",
                                F_IK[ik],I,K,ik,nK,
                                P_JL[jl],J,L,jl,nL);
                      F_contrib(I,J,K,L,i,j,k,l,JK_factor,val,"K ",
                                F_KJ[kj],K,J,kj,nK,
                                P_IL[il],I,L,il,nL);
                      F_contrib(I,J,K,L,i,j,k,l,IL_factor,val,"K ",
                                F_IL[il],I,L,il,nL,
                                P_KJ[kj],K,J,kj,nK);
                      F_contrib(I,J,K,L,i,j,k,l,JL_factor,val,"K ",
                                F_JL[jl],J,L,jl,nL,
                                P_IK[ik],I,K,ik,nK);
                    }
                }
            }
        }
    }
  else {
      double *F_KJ = kmat_block(0, K, J);
      double *F_LJ = kmat_block(0, L, J);
      const double * restrictxx P_KJ = pmat_block(0, K, J);
      const double * restrictxx P_LJ = pmat_block(0, L, J);
      for (int i=0,ijkl=0,il_begin=0,ik_begin=0; i<nI;
           i++,il_begin+=nL,ik_begin+=nK) {
          for (int j=0; j<nJ; j++) {
              for (int k=0,kj=j,ik=ik_begin; k<nK; k++,ik++,kj+=nJ) {
                  for (int l=0, il=il_begin, lj=j;
                       l<nL;
                       l++, ijkl++, il++, lj+=nJ) {
                      double val = buf[ijkl];
                      F_contrib(I,J,K,L,i,j,k,l,IK_factor,val,"K ",
                                F_IK[ik],I,K,ik,nK,
                                P_LJ[lj],L,J,lj,nJ);
                      F_contrib(I,J,K,L,i,j,k,l,JK_factor,val,"K ",
                                F_KJ[kj],K,J,kj,nJ,
                                P_IL[il],I,L,il,nL);
                      F_contrib(I,J,K,L,i,j,k,l,IL_factor,val,"K ",
                                F_IL[il],I,L,il,nL,
                                P_KJ[kj],K,J,kj,nJ);
                      F_contrib(I,J,K,L,i,j,k,l,JL_factor,val,"K ",
                                F_LJ[lj],L,J,lj,nJ,
                                P_IK[ik],I,K,ik,nK);
                    }
                }
            }
        }
    }
}

void
CLHFContribution::contrib_p12_p13p24_J(double factor,
                                       int I, int J, int K, int L,
                                       int nI, int nJ, int nK, int nL,
                                       const double * restrictxx buf)
{
  double factor2 = 2.0 * factor;

  double *F_IJ = jmat_block(0, I, J);
  double *F_KL = jmat_block(0, K, L);
  const double * restrictxx P_IJ = pmat_block(0, I, J);
  const double * restrictxx P_KL = pmat_block(0, K, L);

  int nKL = nK * nL;

  for (int i=0, ijkl=0, ij=0; i<nI; i++) {
      for (int j=0; j<nJ; j++, ij++) {
          double F_IJ_ij = F_IJ[ij];
          double P_IJ_ij = P_IJ[ij];
          for (int kl=0; kl<nKL; kl++, ijkl++) {
              double val = buf[ijkl];
              F_contrib(I,J,K,L,i,j,kl/nL,kl%nL,factor,val,"J ",
                        F_IJ_ij,I,J,ij,nJ,
                        P_KL[kl],K,L,kl,nL);
              F_contrib(I,J,K,L,i,j,kl/nL,kl%nL,factor2,val,"J ",
                        F_KL[kl],K,L,kl,nL,
                        P_IJ_ij,I,J,ij,nJ);
            }
          F_IJ[ij] = F_IJ_ij;
        }
    }
}

void
CLHFContribution::contrib_p12_p13p24_K(double factor,
                                       int I, int J, int K, int L,
                                       int nI, int nJ, int nK, int nL,
                                       const double * restrictxx buf)
{
  double *F_IK = kmat_block(0, I, K);
//   double *F_IL = kmat_block(0, I, L);
//   const double * restrictxx P_IK = pmat_block(0, I, K);
  const double * restrictxx P_IL = pmat_block(0, I, L);

  // We normally leave out the transposed Fock matrix element and apply the
  // standard contribution of -0.5.  However, if both indices are the same,
  // then we need to multiply by two.
  double IK_factor = -0.5*(I==K?2:1) * factor;
  double JK_factor = -0.5*(J==K?2:1) * factor;
//   double IL_factor = -0.5*(I==L?2:1);
//   double JL_factor = -0.5*(J==L?2:1);

  if (J >= K) {
      double *F_JK = kmat_block(0, J, K);
//       double *F_JL = kmat_block(0, J, L);
//       const double * restrictxx P_JK = pmat_block(0, J, K);
      const double * restrictxx P_JL = pmat_block(0, J, L);
      for (int i=0,ijkl=0,il_begin=0,ik_begin=0; i<nI;
           i++,il_begin+=nL,ik_begin+=nK) {
          for (int j=0, jk=0, jl_begin=0; j<nJ; j++, jl_begin += nL) {
              for (int k=0,ik=ik_begin; k<nK; k++,ik++,jk++) {
                  for (int l=0, il=il_begin, jl=jl_begin;
                       l<nL;
                       l++, ijkl++, il++, jl++) {
                      double val = buf[ijkl];
                      F_contrib(I,J,K,L,i,j,k,l,IK_factor,val,"K ",
                                F_IK[ik],I,K,ik,nK,
                                P_JL[jl],J,L,jl,nL);
                      F_contrib(I,J,K,L,i,j,k,l,JK_factor,val,"K ",
                                F_JK[jk],J,K,jk,nK,
                                P_IL[il],I,L,il,nL);
//                       F_contrib(I,J,K,L,i,j,k,l,IL_factor,val,"K ",
//                                 F_IL[il],I,L,il,nL,
//                                 P_JK[jk],J,K,jk,nK);
//                       F_contrib(I,J,K,L,i,j,k,l,JL_factor,val,"K ",
//                                 F_JL[jl],J,L,jl,nL,
//                                 P_IK[ik],I,K,ik,nK);
                    }
                }
            }
        }
    }
  else if (J >= L) {
      double *F_KJ = kmat_block(0, K, J);
//       double *F_JL = kmat_block(0, J, L);
//       const double * restrictxx P_KJ = pmat_block(0, K, J);
      const double * restrictxx P_JL = pmat_block(0, J, L);
      for (int i=0,ijkl=0,il_begin=0,ik_begin=0; i<nI;
           i++,il_begin+=nL,ik_begin+=nK) {
          for (int j=0, jl_begin=0; j<nJ; j++, jl_begin += nL) {
              int kj = j;
              for (int k=0,kj=j,ik=ik_begin; k<nK; k++,ik++,kj+=nJ) {
                  for (int l=0, il=il_begin, jl=jl_begin;
                       l<nL;
                       l++, ijkl++, il++, jl++) {
                      double val = buf[ijkl];
                      F_contrib(I,J,K,L,i,j,k,l,IK_factor,val,"K ",
                                F_IK[ik],I,K,ik,nK,
                                P_JL[jl],J,L,jl,nL);
                      F_contrib(I,J,K,L,i,j,k,l,JK_factor,val,"K ",
                                F_KJ[kj],K,J,kj,nK,
                                P_IL[il],I,L,il,nL);
//                       F_contrib(I,J,K,L,i,j,k,l,IL_factor,val,"K ",
//                                 F_IL[il],I,L,il,nL,
//                                 P_KJ[kj],K,J,kj,nK);
//                       F_contrib(I,J,K,L,i,j,k,l,JL_factor,val,"K ",
//                                 F_JL[jl],J,L,jl,nL,
//                                 P_IK[ik],I,K,ik,nK);
                    }
                }
            }
        }
    }
  else {
      double *F_KJ = kmat_block(0, K, J);
//       double *F_LJ = kmat_block(0, L, J);
//       const double * restrictxx P_KJ = pmat_block(0, K, J);
      const double * restrictxx P_LJ = pmat_block(0, L, J);
      for (int i=0,ijkl=0,il_begin=0,ik_begin=0; i<nI;
           i++,il_begin+=nL,ik_begin+=nK) {
          for (int j=0; j<nJ; j++) {
              for (int k=0,kj=j,ik=ik_begin; k<nK; k++,ik++,kj+=nJ) {
                  for (int l=0, il=il_begin, lj=j;
                       l<nL;
                       l++, ijkl++, il++, lj+=nJ) {
                      double val = buf[ijkl];
                      F_contrib(I,J,K,L,i,j,k,l,IK_factor,val,"K ",
                                F_IK[ik],I,K,ik,nK,
                                P_LJ[lj],L,J,lj,nJ);
                      F_contrib(I,J,K,L,i,j,k,l,JK_factor,val,"K ",
                                F_KJ[kj],K,J,kj,nJ,
                                P_IL[il],I,L,il,nL);
//                       F_contrib(I,J,K,L,i,j,k,l,IL_factor,val,"K ",
//                                 F_IL[il],I,L,il,nL,
//                                 P_KJ[kj],K,J,kj,nJ);
//                       F_contrib(I,J,K,L,i,j,k,l,JL_factor,val,"K ",
//                                 F_LJ[lj],L,J,lj,nJ,
//                                 P_IK[ik],I,K,ik,nK);
                    }
                }
            }
        }
    }
}

void
CLHFContribution::contrib_p34_p13p24_J(double factor,
                                       int I, int J, int K, int L,
                                       int nI, int nJ, int nK, int nL,
                                       const double * restrictxx buf)
{
  double factor2 = 2.0 * factor;

  double *F_IJ = jmat_block(0, I, J);
  double *F_KL = jmat_block(0, K, L);
  const double * restrictxx P_IJ = pmat_block(0, I, J);
  const double * restrictxx P_KL = pmat_block(0, K, L);

  int nKL = nK * nL;

  for (int i=0, ijkl=0, ij=0; i<nI; i++) {
      for (int j=0; j<nJ; j++, ij++) {
          double F_IJ_ij = F_IJ[ij];
          double P_IJ_ij = P_IJ[ij];
          for (int kl=0; kl<nKL; kl++, ijkl++) {
              double val = buf[ijkl];
              F_contrib(I,J,K,L,i,j,kl/nL,kl%nL,factor2,val,"J ",
                        F_IJ_ij,I,J,ij,nJ,
                        P_KL[kl],K,L,kl,nL);
              F_contrib(I,J,K,L,i,j,kl/nL,kl%nL,factor,val,"J ",
                        F_KL[kl],K,L,kl,nL,
                        P_IJ_ij,I,J,ij,nJ);
            }
          F_IJ[ij] = F_IJ_ij;
        }
    }
}

void
CLHFContribution::contrib_p34_p13p24_K(double factor,
                                       int I, int J, int K, int L,
                                       int nI, int nJ, int nK, int nL,
                                       const double * restrictxx buf)
{
  double *F_IK = kmat_block(0, I, K);
  double *F_IL = kmat_block(0, I, L);
  const double * restrictxx P_IK = pmat_block(0, I, K);
  const double * restrictxx P_IL = pmat_block(0, I, L);

  // We normally leave out the transposed Fock matrix element and apply the
  // standard contribution of -0.5.  However, if both indices are the same,
  // then we need to multiply by two.
  double IK_factor = -0.5*(I==K?2:1) * factor;
//   double JK_factor = -0.5*(J==K?2:1);
  double IL_factor = -0.5*(I==L?2:1) * factor;
//   double JL_factor = -0.5*(J==L?2:1);

  if (J >= K) {
//       double *F_JK = kmat_block(0, J, K);
//       double *F_JL = kmat_block(0, J, L);
      const double * restrictxx P_JK = pmat_block(0, J, K);
      const double * restrictxx P_JL = pmat_block(0, J, L);
      for (int i=0,ijkl=0,il_begin=0,ik_begin=0; i<nI;
           i++,il_begin+=nL,ik_begin+=nK) {
          for (int j=0, jk=0, jl_begin=0; j<nJ; j++, jl_begin += nL) {
              for (int k=0,ik=ik_begin; k<nK; k++,ik++,jk++) {
                  for (int l=0, il=il_begin, jl=jl_begin;
                       l<nL;
                       l++, ijkl++, il++, jl++) {
                      double val = buf[ijkl];
                      F_contrib(I,J,K,L,i,j,k,l,IK_factor,val,"K ",
                                F_IK[ik],I,K,ik,nK,
                                P_JL[jl],J,L,jl,nL);
//                       F_contrib(I,J,K,L,i,j,k,l,JK_factor,val,"K ",
//                                 F_JK[jk],J,K,jk,nK,
//                                 P_IL[il],I,L,il,nL);
                      F_contrib(I,J,K,L,i,j,k,l,IL_factor,val,"K ",
                                F_IL[il],I,L,il,nL,
                                P_JK[jk],J,K,jk,nK);
//                       F_contrib(I,J,K,L,i,j,k,l,JL_factor,val,"K ",
//                                 F_JL[jl],J,L,jl,nL,
//                                 P_IK[ik],I,K,ik,nK);
                    }
                }
            }
        }
    }
  else if (J >= L) {
//       double *F_KJ = kmat_block(0, K, J);
//       double *F_JL = kmat_block(0, J, L);
      const double * restrictxx P_KJ = pmat_block(0, K, J);
      const double * restrictxx P_JL = pmat_block(0, J, L);
      for (int i=0,ijkl=0,il_begin=0,ik_begin=0; i<nI;
           i++,il_begin+=nL,ik_begin+=nK) {
          for (int j=0, jl_begin=0; j<nJ; j++, jl_begin += nL) {
              int kj = j;
              for (int k=0,kj=j,ik=ik_begin; k<nK; k++,ik++,kj+=nJ) {
                  for (int l=0, il=il_begin, jl=jl_begin;
                       l<nL;
                       l++, ijkl++, il++, jl++) {
                      double val = buf[ijkl];
                      F_contrib(I,J,K,L,i,j,k,l,IK_factor,val,"K ",
                                F_IK[ik],I,K,ik,nK,
                                P_JL[jl],J,L,jl,nL);
//                       F_contrib(I,J,K,L,i,j,k,l,JK_factor,val,"K ",
//                                 F_KJ[kj],K,J,kj,nK,
//                                 P_IL[il],I,L,il,nL);
                      F_contrib(I,J,K,L,i,j,k,l,IL_factor,val,"K ",
                                F_IL[il],I,L,il,nL,
                                P_KJ[kj],K,J,kj,nK);
//                       F_contrib(I,J,K,L,i,j,k,l,JL_factor,val,"K ",
//                                 F_JL[jl],J,L,jl,nL,
//                                 P_IK[ik],I,K,ik,nK);
                    }
                }
            }
        }
    }
  else {
//       double *F_KJ = kmat_block(0, K, J);
//       double *F_LJ = kmat_block(0, L, J);
      const double * restrictxx P_KJ = pmat_block(0, K, J);
      const double * restrictxx P_LJ = pmat_block(0, L, J);
      for (int i=0,ijkl=0,il_begin=0,ik_begin=0; i<nI;
           i++,il_begin+=nL,ik_begin+=nK) {
          for (int j=0; j<nJ; j++) {
              for (int k=0,kj=j,ik=ik_begin; k<nK; k++,ik++,kj+=nJ) {
                  for (int l=0, il=il_begin, lj=j;
                       l<nL;
                       l++, ijkl++, il++, lj+=nJ) {
                      double val = buf[ijkl];
                      F_contrib(I,J,K,L,i,j,k,l,IK_factor,val,"K ",
                                F_IK[ik],I,K,ik,nK,
                                P_LJ[lj],L,J,lj,nJ);
//                       F_contrib(I,J,K,L,i,j,k,l,JK_factor,val,"K ",
//                                 F_KJ[kj],K,J,kj,nJ,
//                                 P_IL[il],I,L,il,nL);
                      F_contrib(I,J,K,L,i,j,k,l,IL_factor,val,"K ",
                                F_IL[il],I,L,il,nL,
                                P_KJ[kj],K,J,kj,nJ);
//                       F_contrib(I,J,K,L,i,j,k,l,JL_factor,val,"K ",
//                                 F_LJ[lj],L,J,lj,nJ,
//                                 P_IK[ik],I,K,ik,nK);
                    }
                }
            }
        }
    }
}

void
CLHFContribution::contrib_p12_p34_J(double factor,
                                    int I, int J, int K, int L,
                                    int nI, int nJ, int nK, int nL,
                                    const double * restrictxx buf)
{
  double factor2 = 2.0 * factor;

  double *F_IJ = jmat_block(0, I, J);
  const double * restrictxx P_KL = pmat_block(0, K, L);

  int nKL = nK * nL;

  for (int i=0, ijkl=0, ij=0; i<nI; i++) {
      for (int j=0; j<nJ; j++, ij++) {
          double F_IJ_ij = F_IJ[ij];
          for (int kl=0; kl<nKL; kl++, ijkl++) {
              double val = buf[ijkl];
              F_contrib(I,J,K,L,i,j,kl/nL,kl%nL,factor2,val,"J ",
                        F_IJ_ij,I,J,ij,nJ,
                        P_KL[kl],K,L,kl,nL);
            }
          F_IJ[ij] = F_IJ_ij;
        }
    }
}

void
CLHFContribution::contrib_p12_p34_K(double factor,
                                    int I, int J, int K, int L,
                                    int nI, int nJ, int nK, int nL,
                                    const double * restrictxx buf)
{
  double *F_IK = kmat_block(0, I, K);
  double *F_IL = kmat_block(0, I, L);
  const double * restrictxx P_IK = pmat_block(0, I, K);
  const double * restrictxx P_IL = pmat_block(0, I, L);


  double IK_factor = -0.5*factor; // IK are always canonical
  double IL_factor = -0.5*factor; // IK are always canonical
  double JK_factor = -0.5*(J>=K?1:0)*factor;
  double JL_factor = -0.5*factor; // JL are always canonical

  if (J >= K) {
      double *F_JK = kmat_block(0, J, K);
      double *F_JL = kmat_block(0, J, L);
      const double * restrictxx P_JK = pmat_block(0, J, K);
      const double * restrictxx P_JL = pmat_block(0, J, L);
      for (int i=0,ijkl=0,il_begin=0,ik_begin=0; i<nI;
           i++,il_begin+=nL,ik_begin+=nK) {
          for (int j=0, jk=0, jl_begin=0; j<nJ; j++, jl_begin += nL) {
              for (int k=0,ik=ik_begin; k<nK; k++,ik++,jk++) {
                  for (int l=0, il=il_begin, jl=jl_begin;
                       l<nL;
                       l++, ijkl++, il++, jl++) {
                      double val = buf[ijkl];
                      F_contrib(I,J,K,L,i,j,k,l,IK_factor,val,"K ",
                                F_IK[ik],I,K,ik,nK,
                                P_JL[jl],J,L,jl,nL);
                      F_contrib(I,J,K,L,i,j,k,l,JK_factor,val,"K ",
                                F_JK[jk],J,K,jk,nK,
                                P_IL[il],I,L,il,nL);
                      F_contrib(I,J,K,L,i,j,k,l,IL_factor,val,"K ",
                                F_IL[il],I,L,il,nL,
                                P_JK[jk],J,K,jk,nK);
                      F_contrib(I,J,K,L,i,j,k,l,JL_factor,val,"K ",
                                F_JL[jl],J,L,jl,nL,
                                P_IK[ik],I,K,ik,nK);
                    }
                }
            }
        }
    }
  else if (J >= L) {
      double *F_KJ = kmat_block(0, K, J);
      double *F_JL = kmat_block(0, J, L);
      const double * restrictxx P_KJ = pmat_block(0, K, J);
      const double * restrictxx P_JL = pmat_block(0, J, L);
      for (int i=0,ijkl=0,il_begin=0,ik_begin=0; i<nI;
           i++,il_begin+=nL,ik_begin+=nK) {
          for (int j=0, jl_begin=0; j<nJ; j++, jl_begin += nL) {
              int kj = j;
              for (int k=0,kj=j,ik=ik_begin; k<nK; k++,ik++,kj+=nJ) {
                  for (int l=0, il=il_begin, jl=jl_begin;
                       l<nL;
                       l++, ijkl++, il++, jl++) {
                      double val = buf[ijkl];
                      F_contrib(I,J,K,L,i,j,k,l,IK_factor,val,"K ",
                                F_IK[ik],I,K,ik,nK,
                                P_JL[jl],J,L,jl,nL);
                      F_contrib(I,J,K,L,i,j,k,l,JK_factor,val,"K ",
                                F_KJ[kj],K,J,kj,nK,
                                P_IL[il],I,L,il,nL);
                      F_contrib(I,J,K,L,i,j,k,l,IL_factor,val,"K ",
                                F_IL[il],I,L,il,nL,
                                P_KJ[kj],K,J,kj,nK);
                      F_contrib(I,J,K,L,i,j,k,l,JL_factor,val,"K ",
                                F_JL[jl],J,L,jl,nL,
                                P_IK[ik],I,K,ik,nK);
                    }
                }
            }
        }
    }
  else {
      double *F_KJ = kmat_block(0, K, J);
      double *F_LJ = kmat_block(0, L, J);
      const double * restrictxx P_KJ = pmat_block(0, K, J);
      const double * restrictxx P_LJ = pmat_block(0, L, J);
      for (int i=0,ijkl=0,il_begin=0,ik_begin=0; i<nI;
           i++,il_begin+=nL,ik_begin+=nK) {
          for (int j=0; j<nJ; j++) {
              for (int k=0,kj=j,ik=ik_begin; k<nK; k++,ik++,kj+=nJ) {
                  for (int l=0, il=il_begin, lj=j;
                       l<nL;
                       l++, ijkl++, il++, lj+=nJ) {
                      double val = buf[ijkl];
                      F_contrib(I,J,K,L,i,j,k,l,IK_factor,val,"K ",
                                F_IK[ik],I,K,ik,nK,
                                P_LJ[lj],L,J,lj,nJ);
                      F_contrib(I,J,K,L,i,j,k,l,JK_factor,val,"K ",
                                F_KJ[kj],K,J,kj,nJ,
                                P_IL[il],I,L,il,nL);
                      F_contrib(I,J,K,L,i,j,k,l,IL_factor,val,"K ",
                                F_IL[il],I,L,il,nL,
                                P_KJ[kj],K,J,kj,nJ);
                      F_contrib(I,J,K,L,i,j,k,l,JL_factor,val,"K ",
                                F_LJ[lj],L,J,lj,nJ,
                                P_IK[ik],I,K,ik,nK);
                    }
                }
            }
        }
    }
}

void
CLHFContribution::contrib_p12_J(double factor,
                                int I, int J, int K, int L,
                                int nI, int nJ, int nK, int nL,
                                const double * restrictxx buf)
{
  contrib_e_J(factor,I,J,K,L,nI,nJ,nK,nL,buf);
}

void
CLHFContribution::contrib_p12_K(double factor,
                                int I, int J, int K, int L,
                                int nI, int nJ, int nK, int nL,
                                const double * restrictxx buf)
{
  double *F_IK = kmat_block(0, I, K);
  const double * restrictxx P_IL = pmat_block(0, I, L);

  double IK_factor = -0.5*factor;
  double JK_factor = -0.5*factor;

  double *F_JK = kmat_block(0, J, K);
  const double * restrictxx P_JL = pmat_block(0, J, L);
  for (int i=0,ijkl=0,il_begin=0,ik_begin=0; i<nI;
       i++,il_begin+=nL,ik_begin+=nK) {
      for (int j=0, jk=0, jl_begin=0; j<nJ; j++, jl_begin += nL) {
          for (int k=0,ik=ik_begin; k<nK; k++,ik++,jk++) {
              for (int l=0, il=il_begin, jl=jl_begin;
                   l<nL;
                   l++, ijkl++, il++, jl++) {
                  double val = buf[ijkl];
                  F_contrib(I,J,K,L,i,j,k,l,IK_factor,val,"K ",
                            F_IK[ik],I,K,ik,nK,
                            P_JL[jl],J,L,jl,nL);
                  F_contrib(I,J,K,L,i,j,k,l,JK_factor,val,"K ",
                            F_JK[jk],J,K,jk,nK,
                            P_IL[il],I,L,il,nL);
                }
            }
        }
    }
}

void
CLHFContribution::contrib_p34_J(double factor,
                                int I, int J, int K, int L,
                                int nI, int nJ, int nK, int nL,
                                const double * restrictxx buf)
{
  contrib_p12_p34_J(factor,I,J,K,L,nI,nJ,nK,nL,buf);
}

void
CLHFContribution::contrib_p34_K(double factor,
                                int I, int J, int K, int L,
                                int nI, int nJ, int nK, int nL,
                                const double * restrictxx buf)
{
  double *F_IK = kmat_block(0, I, K);
  double *F_IL = kmat_block(0, I, L);

  double IK_factor = -0.5*factor;
  double IL_factor = -0.5*factor;

  double *F_JK = kmat_block(0, J, K);
  double *F_JL = kmat_block(0, J, L);
  const double * restrictxx P_JK = pmat_block(0, J, K);
  const double * restrictxx P_JL = pmat_block(0, J, L);
  for (int i=0,ijkl=0,il_begin=0,ik_begin=0; i<nI;
       i++,il_begin+=nL,ik_begin+=nK) {
      for (int j=0, jk=0, jl_begin=0; j<nJ; j++, jl_begin += nL) {
          for (int k=0,ik=ik_begin; k<nK; k++,ik++,jk++) {
              for (int l=0, il=il_begin, jl=jl_begin;
                   l<nL;
                   l++, ijkl++, il++, jl++) {
                  double val = buf[ijkl];
                  F_contrib(I,J,K,L,i,j,k,l,IK_factor,val,"K ",
                            F_IK[ik],I,K,ik,nK,
                            P_JL[jl],J,L,jl,nL);
                  F_contrib(I,J,K,L,i,j,k,l,IL_factor,val,"K ",
                            F_IL[il],I,L,il,nL,
                            P_JK[jk],J,K,jk,nK);
                }
            }
        }
    }
}

void
CLHFContribution::contrib_p13p24_J(double factor,
                                   int I, int J, int K, int L,
                                   int nI, int nJ, int nK, int nL,
                                   const double * restrictxx buf)
{
  double *F_IJ = jmat_block(0, I, J);
  double *F_KL = jmat_block(0, K, L);
  const double * restrictxx P_IJ = pmat_block(0, I, J);
  const double * restrictxx P_KL = pmat_block(0, K, L);

  int nKL = nK * nL;

  for (int i=0, ijkl=0, ij=0; i<nI; i++) {
      for (int j=0; j<nJ; j++, ij++) {
          double F_IJ_ij = F_IJ[ij];
          double P_IJ_ij = P_IJ[ij];
          for (int kl=0; kl<nKL; kl++, ijkl++) {
              double val = buf[ijkl];
              F_contrib(I,J,K,L,i,j,kl/nL,kl%nL,factor,val,"J ",
                        F_IJ_ij,I,J,ij,nJ,
                        P_KL[kl],K,L,kl,nL);
              F_contrib(I,J,K,L,i,j,kl/nL,kl%nL,factor,val,"J ",
                        F_KL[kl],K,L,kl,nL,
                        P_IJ_ij,I,J,ij,nJ);
            }
          F_IJ[ij] = F_IJ_ij;
        }
    }
}

void
CLHFContribution::contrib_p13p24_K(double factor,
                                   int I, int J, int K, int L,
                                   int nI, int nJ, int nK, int nL,
                                   const double * restrictxx buf)
{
  double *F_IK = kmat_block(0, I, K);
//   double *F_IL = kmat_block(0, I, L);
//   const double * restrictxx P_IK = pmat_block(0, I, K);
//   const double * restrictxx P_IL = pmat_block(0, I, L);

  // We normally leave out the transposed Fock matrix element and apply the
  // standard contribution of -0.5.  However, if both indices are the same,
  // then we need to multiply by two.
  double IK_factor = -0.5*(I==K?2:1)*factor;
//   double JK_factor = -0.5*(J==K?2:1);
//   double IL_factor = -0.5*(I==L?2:1);
//   double JL_factor = -0.5*(J==L?2:1);

  if (J >= K) {
//       double *F_JK = kmat_block(0, J, K);
//       double *F_JL = kmat_block(0, J, L);
//       const double * restrictxx P_JK = pmat_block(0, J, K);
      const double * restrictxx P_JL = pmat_block(0, J, L);
      for (int i=0,ijkl=0,il_begin=0,ik_begin=0; i<nI;
           i++,il_begin+=nL,ik_begin+=nK) {
          for (int j=0, jk=0, jl_begin=0; j<nJ; j++, jl_begin += nL) {
              for (int k=0,ik=ik_begin; k<nK; k++,ik++,jk++) {
                  for (int l=0, il=il_begin, jl=jl_begin;
                       l<nL;
                       l++, ijkl++, il++, jl++) {
                      double val = buf[ijkl];
                      F_contrib(I,J,K,L,i,j,k,l,IK_factor,val,"K ",
                                F_IK[ik],I,K,ik,nK,
                                P_JL[jl],J,L,jl,nL);
//                       F_contrib(I,J,K,L,i,j,k,l,JK_factor,val,"K ",
//                                 F_JK[jk],J,K,jk,nK,
//                                 P_IL[il],I,L,il,nL);
//                       F_contrib(I,J,K,L,i,j,k,l,IL_factor,val,"K ",
//                                 F_IL[il],I,L,il,nL,
//                                 P_JK[jk],J,K,jk,nK);
//                       F_contrib(I,J,K,L,i,j,k,l,JL_factor,val,"K ",
//                                 F_JL[jl],J,L,jl,nL,
//                                 P_IK[ik],I,K,ik,nK);
                    }
                }
            }
        }
    }
  else if (J >= L) {
//       double *F_KJ = kmat_block(0, K, J);
//       double *F_JL = kmat_block(0, J, L);
//       const double * restrictxx P_KJ = pmat_block(0, K, J);
      const double * restrictxx P_JL = pmat_block(0, J, L);
      for (int i=0,ijkl=0,il_begin=0,ik_begin=0; i<nI;
           i++,il_begin+=nL,ik_begin+=nK) {
          for (int j=0, jl_begin=0; j<nJ; j++, jl_begin += nL) {
              int kj = j;
              for (int k=0,kj=j,ik=ik_begin; k<nK; k++,ik++,kj+=nJ) {
                  for (int l=0, il=il_begin, jl=jl_begin;
                       l<nL;
                       l++, ijkl++, il++, jl++) {
                      double val = buf[ijkl];
                      F_contrib(I,J,K,L,i,j,k,l,IK_factor,val,"K ",
                                F_IK[ik],I,K,ik,nK,
                                P_JL[jl],J,L,jl,nL);
//                       F_contrib(I,J,K,L,i,j,k,l,JK_factor,val,"K ",
//                                 F_KJ[kj],K,J,kj,nK,
//                                 P_IL[il],I,L,il,nL);
//                       F_contrib(I,J,K,L,i,j,k,l,IL_factor,val,"K ",
//                                 F_IL[il],I,L,il,nL,
//                                 P_KJ[kj],K,J,kj,nK);
//                       F_contrib(I,J,K,L,i,j,k,l,JL_factor,val,"K ",
//                                 F_JL[jl],J,L,jl,nL,
//                                 P_IK[ik],I,K,ik,nK);
                    }
                }
            }
        }
    }
  else {
//       double *F_KJ = kmat_block(0, K, J);
//       double *F_LJ = kmat_block(0, L, J);
//       const double * restrictxx P_KJ = pmat_block(0, K, J);
      const double * restrictxx P_LJ = pmat_block(0, L, J);
      for (int i=0,ijkl=0,il_begin=0,ik_begin=0; i<nI;
           i++,il_begin+=nL,ik_begin+=nK) {
          for (int j=0; j<nJ; j++) {
              for (int k=0,kj=j,ik=ik_begin; k<nK; k++,ik++,kj+=nJ) {
                  for (int l=0, il=il_begin, lj=j;
                       l<nL;
                       l++, ijkl++, il++, lj+=nJ) {
                      double val = buf[ijkl];
                      F_contrib(I,J,K,L,i,j,k,l,IK_factor,val,"K ",
                                F_IK[ik],I,K,ik,nK,
                                P_LJ[lj],L,J,lj,nJ);
//                       F_contrib(I,J,K,L,i,j,k,l,JK_factor,val,"K ",
//                                 F_KJ[kj],K,J,kj,nJ,
//                                 P_IL[il],I,L,il,nL);
//                       F_contrib(I,J,K,L,i,j,k,l,IL_factor,val,"K ",
//                                 F_IL[il],I,L,il,nL,
//                                 P_KJ[kj],K,J,kj,nJ);
//                       F_contrib(I,J,K,L,i,j,k,l,JL_factor,val,"K ",
//                                 F_LJ[lj],L,J,lj,nJ,
//                                 P_IK[ik],I,K,ik,nK);
                    }
                }
            }
        }
    }
}

};

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
