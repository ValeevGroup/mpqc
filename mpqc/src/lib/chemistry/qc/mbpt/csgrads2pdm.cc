//
// csgrads2pdm.cc
// based on csgrad.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Ida Nielsen <ida@kemi.aau.dk>
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include <util/misc/formio.h>
#include <util/group/message.h>
#include <chemistry/qc/mbpt/csgrads2pdm.h>

using namespace sc;

CSGradS2PDM::CSGradS2PDM(int mythread_a, int nthread_a,
                         int me_a, int nproc_a,
                         const Ref<ThreadLock> &lock_a,
                         const Ref<GaussianBasisSet> &basis_a,
                         const Ref<TwoBodyDerivInt> &tbintder_a,
                         const double *PHF_a, const double *P2AO_a,
                         int tol_a, int debug_a, int dynamic_a)
{
  mythread = mythread_a;
  nthread = nthread_a;
  me = me_a;
  nproc = nproc_a;
  lock = lock_a;
  basis = basis_a;
  tbintder = tbintder_a;
  PHF = PHF_a;
  P2AO = P2AO_a;
  tol = tol_a;
  debug = debug_a;
  dynamic = dynamic_a;

  int natom = basis->molecule()->natom();
  ginter = new double*[natom];
  ginter[0] = new double[natom*3];
  hf_ginter = new double*[natom];
  hf_ginter[0] = new double[natom*3];
  for (int i=0; i<natom; i++) {
    ginter[i] = &ginter[0][i*3];
    hf_ginter[i] = &hf_ginter[0][i*3];
    for (int j=0; j<3; j++) {
      ginter[i][j] = 0.0;
      hf_ginter[i][j] = 0.0;
      }
    }
}

CSGradS2PDM::~CSGradS2PDM()
{
  delete[] ginter[0];
  delete[] ginter;
  delete[] hf_ginter[0];
  delete[] hf_ginter;
}

void
CSGradS2PDM::accum_contrib(double **sum, double **contribs)
{
  int natom = basis->molecule()->natom();
  for (int i=0; i<natom; i++) {
    for (int j=0; j<3; j++) {
      sum[i][j] += contribs[i][j];
      }
    }
}

void
CSGradS2PDM::accum_mp2_contrib(double **ginter_a)
{
  accum_contrib(ginter_a, ginter);
}

void
CSGradS2PDM::accum_hf_contrib(double **hf_ginter_a)
{
  accum_contrib(hf_ginter_a, hf_ginter);
}

//////////////////////////////////////////////////////////////
// Compute (in the AO basis) the contribution to the gradient 
// from the separable part of the two particle density matrix
//////////////////////////////////////////////////////////////
void
CSGradS2PDM::run()
{               

  int P, Q, R, S;
  int QP, SR;
  int p, q, r;
  int np, nq, nr, ns;
  int p_offset, q_offset, r_offset, s_offset;
  int bf1, bf2, bf3, bf4;
  int xyz;
  int derset;
  int nshell = basis->nshell();
  int nbasis = basis->nbasis();

  double *grad_ptr1, *grad_ptr2;
  double *hf_grad_ptr1, *hf_grad_ptr2;
  double tmpval;
  const double *phf_pq, *phf_pr, *phf_ps, *phf_qr, *phf_qs, *phf_rs;
  const double *p2ao_pq, *p2ao_pr, *p2ao_ps, *p2ao_qr, *p2ao_qs, *p2ao_rs;
  double k_QP, k_SR, k_QPSR; // factors needed since we loop over nonredund
                             // shell quartets but do redund integrals within
                             // shell quartets when applicable
  double gamma_factor; // factor multiplying integrals; needed because we
                       // loop over nonredund shell quarters but do redund
                       // integrals within shell quartets when applicable
  double *gammasym_pqrs; // symmetrized sep. 2PDM
  double *gammasym_ptr;
  double *hf_gammasym_pqrs; // HF only versions of gammsym
  double *hf_gammasym_ptr;
  const double *integral_ptr;
  int nfuncmax = basis->max_nfunction_in_shell();
  const double *intderbuf = tbintder->buffer();

  gammasym_pqrs = new double[nfuncmax*nfuncmax*nfuncmax*nfuncmax];
  hf_gammasym_pqrs = new double[nfuncmax*nfuncmax*nfuncmax*nfuncmax];

  DerivCenters der_centers;

  int index = 0;
  int threadindex = 0;

  for (Q=0; Q<nshell; Q++) {
    nq = basis->shell(Q).nfunction();
    q_offset = basis->shell_to_function(Q);

    for (S=0; S<=Q; S++) {
      ns = basis->shell(S).nfunction();
      s_offset = basis->shell_to_function(S);

      for (R=0; R<=S; R++) {
        nr = basis->shell(R).nfunction();
        r_offset = basis->shell_to_function(R);
        k_SR = (R == S ? 0.5 : 1.0);
        SR = S*(S+1)/2 + R;

        for (P=0; P<=(S==Q ? R:Q); P++) {
          // If integral derivative is 0, skip to next P
          if (tbintder->log2_shell_bound(P,Q,R,S) < tol) continue;

          if (index++%nproc == me && threadindex++%nthread == mythread) {
            np = basis->shell(P).nfunction();
            p_offset = basis->shell_to_function(P);
            k_QP = (P == Q ? 0.5 : 1.0);
            QP = Q*(Q+1)/2 + P;
            k_QPSR = (QP == SR ? 0.5 : 1.0);
            gamma_factor = k_QP*k_SR*k_QPSR;

            // Evaluate derivative integrals
            tbintder->compute_shell(P,Q,R,S,der_centers);

            //////////////////////////////
            // Symmetrize sep. 2PDM
            //////////////////////////////
            gammasym_ptr = gammasym_pqrs;
            hf_gammasym_ptr = hf_gammasym_pqrs;
            for (bf1=0; bf1<np; bf1++) {
              p = p_offset + bf1;
              phf_pq =  &PHF [p*nbasis + q_offset];
              p2ao_pq = &P2AO[p*nbasis + q_offset];

              for (bf2=0; bf2<nq; bf2++) {
                q = q_offset + bf2;
                phf_pr =  &PHF [p*nbasis + r_offset];
                p2ao_pr = &P2AO[p*nbasis + r_offset];
                phf_qr =  &PHF [q*nbasis + r_offset];
                p2ao_qr = &P2AO[q*nbasis + r_offset];

                for (bf3=0; bf3<nr; bf3++) {
                  r = r_offset + bf3;
                  phf_ps =  &PHF [p*nbasis + s_offset];
                  phf_qs =  &PHF [q*nbasis + s_offset];
                  phf_rs =  &PHF [r*nbasis + s_offset];
                  p2ao_ps = &P2AO[p*nbasis + s_offset];
                  p2ao_qs = &P2AO[q*nbasis + s_offset];
                  p2ao_rs = &P2AO[r*nbasis + s_offset];

                  for (bf4=0; bf4<ns; bf4++) {

                    *gammasym_ptr++ = gamma_factor*(
                                          4**phf_pq*(*phf_rs + *p2ao_rs)
                                        + 4**phf_rs**p2ao_pq
                                        - *phf_qs*(*phf_pr + *p2ao_pr)
                                        - *phf_qr*(*phf_ps + *p2ao_ps)
                                        - *phf_ps**p2ao_qr
                                        - *phf_pr**p2ao_qs);

                    *hf_gammasym_ptr++ = gamma_factor*(
                                          4**phf_pq*(*phf_rs)
                                        - *phf_qs*(*phf_pr)
                                        - *phf_qr*(*phf_ps));

                    phf_ps++;
                    phf_qs++;
                    phf_rs++;
                    p2ao_ps++;
                    p2ao_qs++;
                    p2ao_rs++;
                    } // exit bf4 loop
                  phf_pr++;
                  p2ao_pr++;
                  phf_qr++;
                  p2ao_qr++;
                  }   // exit bf3 loop
                phf_pq++;
                p2ao_pq++;
                }     // exit bf2 loop
              }       // exit bf1 loop

            ///////////////////////////////////////////////////////////
            // Contract symmetrized sep 2PDM with integral derivatives
            ///////////////////////////////////////////////////////////
            integral_ptr = intderbuf;
            for (derset=0; derset<der_centers.n(); derset++) {

              for (xyz=0; xyz<3; xyz++) {
                grad_ptr1 = &ginter[der_centers.atom(derset)][xyz];
                hf_grad_ptr1 = &hf_ginter[der_centers.atom(derset)][xyz];
                if (der_centers.has_omitted_center()) {
                  grad_ptr2 = &ginter[der_centers.omitted_atom()][xyz];
                  hf_grad_ptr2 = &hf_ginter[der_centers.omitted_atom()][xyz];
                  }

                gammasym_ptr = gammasym_pqrs;
                hf_gammasym_ptr = hf_gammasym_pqrs;
                for (bf1=0; bf1<np; bf1++) {
                  for (bf2=0; bf2<nq; bf2++) {
                    for (bf3=0; bf3<nr; bf3++) {
                      for (bf4=0; bf4<ns; bf4++) {
                        double intval = *integral_ptr++;
                        tmpval = intval * *gammasym_ptr++;
                        *grad_ptr1 += tmpval;
                        *grad_ptr2 -= tmpval;
                        tmpval = intval * *hf_gammasym_ptr++;
                        *hf_grad_ptr1 += tmpval;
                        *hf_grad_ptr2 -= tmpval;
                        } // exit bf4 loop
                      }   // exit bf3 loop
                    }     // exit bf2 loop
                  }       // exit bf1 loop
                }         // exit xyz loop
              }           // exit derset loop


            } // exit "if"
          }   // exit P loop
        }     // exit R loop
      }       // exit S loop
    }         // exit Q loop

  delete[] gammasym_pqrs;
  delete[] hf_gammasym_pqrs;

}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
