//
// hsosv1e1.cc
// based on: csgrade12.cc
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

#include <math.h>

#include <util/misc/formio.h>
#include <chemistry/qc/mbpt/bzerofast.h>
#include <chemistry/qc/mbpt/hsosv1e1.h>

using namespace sc;

#define PRINT1Q 0

HSOSV1Erep1Qtr::HSOSV1Erep1Qtr(int mythread_a, int nthread_a,
                               int me_a, int nproc_a,
                               const Ref<ThreadLock> &lock_a,
                               const Ref<GaussianBasisSet> &basis_a,
                               const Ref<TwoBodyInt> &tbint_a,
                               int ni_a, double **scf_vector_a,
                               double tol_a, int debug_a)
{
  mythread = mythread_a;
  nthread = nthread_a;
  me = me_a;
  nproc = nproc_a;
  lock = lock_a;
  basis = basis_a;
  tbint = tbint_a;
  ni = ni_a;
  scf_vector = scf_vector_a;
  tol = tol_a;
  debug = debug_a;

  nbasis = basis->nbasis();
  nfuncmax = basis->max_nfunction_in_shell();
  nshell = basis->nshell();

  aoint_computed_ = 0;

  trans_int1 = new double[nfuncmax*nfuncmax*nbasis*ni];

  timer = new RegionTimer();
}

HSOSV1Erep1Qtr::~HSOSV1Erep1Qtr()
{
  delete[] trans_int1;
}

void
HSOSV1Erep1Qtr::accum_buffer(double *buffer)
{
  int n = nr*ns*nbasis*ni;
  for (int i=0; i<n; i++) {
    buffer[i] += trans_int1[i];
    }
}

void
HSOSV1Erep1Qtr::set_data(int R_a,int nr_a,int S_a,int ns_a,int ni_a,int ioffset_a)
{
  R = R_a;
  nr = nr_a;
  S = S_a;
  ns = ns_a;
  ni = ni_a;
  i_offset = ioffset_a;
}

void
HSOSV1Erep1Qtr::run()
{
  int i;
  int P, Q;
  int bf1,bf2,bf3,bf4;
  int p,q;
  double *c_pi,*c_qi;

  Timer tim(timer);
  tim.enter("bzerofast trans_int1");
  bzerofast(trans_int1,nfuncmax*nfuncmax*nbasis*ni);
  tim.exit("bzerofast trans_int1");

  const double *intbuf = tbint->buffer();

  int shell_index = 0;
  int thindex = 0;

  for (P = 0; P < basis->nshell(); P++) {
    int np = basis->shell(P).nfunction();

    for (Q = 0; Q <= P; Q++) {
      shell_index++;
      if (shell_index%nproc != me) continue; 
      if (thindex++%nthread != mythread) continue;

      if (tbint->log2_shell_bound(P,Q,R,S) < tol) {
        continue;                           /* skip ereps less than tol */
        }

      aoint_computed_++;

      int nq = basis->shell(Q).nfunction();

      tim.enter("erep");
      tbint->compute_shell(P,Q,R,S);
      tim.exit("erep");

      tim.enter("1. quart. tr."); 

      int index = 0;

      for (bf1 = 0; bf1 < np; bf1++) {
        p = basis->shell_to_function(P) + bf1;
 
        for (bf2 = 0; bf2 < nq; bf2++) {
          q = basis->shell_to_function(Q) + bf2;
          if (q > p) {
            /* if q > p: want to skip the loops over bf3-4  */
            /* and larger bf2 values, so increment bf1 by 1 */
            /* ("break") and adjust the value of index      */
            index = (bf1 + 1) * nq * nr * ns;
            break;
            }

          for (bf3 = 0; bf3 < nr; bf3++) {

            for (bf4 = 0; bf4 < ns; bf4++,index++) {
              if (R==S && bf4>bf3) {
                index = ((bf1*nq + bf2)*nr + (bf3+1))*ns;
                break; 
                }

              if (fabs(intbuf[index])>1.0e-15) {
                double pqrs = intbuf[index];

                double *iqrs = &trans_int1[((bf4*nr + bf3)*nbasis + q)*ni];
                double *iprs = &trans_int1[((bf4*nr + bf3)*nbasis + p)*ni];
                    
                if (p == q) pqrs *= 0.5;

                int col_index = i_offset;
                c_pi = &scf_vector[p][col_index];
                c_qi = &scf_vector[q][col_index];

                for (i=ni; i; i--) {
                  *iqrs++ += pqrs * *c_pi++;
                  *iprs++ += pqrs * *c_qi++;
                  }
                }
              }   /* exit bf4 loop */
            }     /* exit bf3 loop */
          }       /* exit bf2 loop */
        }         /* exit bf1 loop */
      tim.exit("1. quart. tr.");
      }           /* exit Q loop */
    }             /* exit P loop */
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
