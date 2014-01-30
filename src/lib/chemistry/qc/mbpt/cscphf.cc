//
// cscphf.cc
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

#include <math.h>
#include <stdlib.h>
#include <iostream>

#include <util/misc/formio.h>
#include <util/keyval/keyval.h>
#include <math/scmat/matrix.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/mbpt/mbpt.h>

using namespace std;
using namespace sc;

static void
compute_alpha(int dim, double **AP, double **alpha,
              double **P, double *eigval, int nocc, int nvir);

//////////////////////////////////////////////////////////////////////////
// Do a direct CPHF calculation in the AO basis; equations are formulated
// in terms of the occ-vir block P2aj of the second order correction (P2) to 
// the MP2 density matrix (cf. Frisch et al., CPL 166, p. 275 (1990)).
//
// CPHF equations: 
// (I-A)P2aj - B = 0 (B(a,j) = L(a,j)/(eigval[a]-eigval[j]))
// A is a matrix (dimension dimP*dimP),
// P2aj and B are vectors (dimension dimP)
//   (P2aj is kept as a RefSCMatrix);
// Only closed-shell cases handled; no orbitals can be frozen
// On exit, P2aj has been computed.

void 
MBPT2::cs_cphf(double **scf_vector,
               double *Laj, double *eigval, RefSCMatrix& P2aj)
{
  double epsilon = cphf_epsilon_; //convergence criterion for P2aj

  int i, j, k, l, a;
  int niter;
  int dimP = nocc*nvir;

  Ref<SCMatrixKit> kit = basis()->matrixkit();

  RefSCDimension nbasis_dim = ao_dimension()->blocks()->subdim(0);
  RefSCDimension nvir_dim(new SCDimension(nvir,1));
  nvir_dim->blocks()->set_subdim(0, new SCDimension(nvir));
  RefSCDimension nocc_dim(new SCDimension(nocc,1));
  nocc_dim->blocks()->set_subdim(0, new SCDimension(nocc));

  RefSCMatrix Cv(nbasis_dim,nvir_dim,kit);
  RefSCMatrix Co(nbasis_dim,nocc_dim,kit);
  RefSCMatrix D_matrix(nbasis_dim,nbasis_dim,kit);
  RefSCMatrix AP_matrix(nvir_dim,nocc_dim,kit); // holds A*P[i-1]
  RefSCMatrix P_matrix(nvir_dim, nocc_dim, kit);
  RefSymmSCMatrix G(nbasis_dim,kit);


  double *projctn = new double[dimP];
  double *P_sum_new = new double[dimP];
  double *P_sum_old = new double[dimP];
  double **AP_matrix_tot; // row is A*P[k]
  double **P_tmp, **alpha_tmp, **AP_matrix_tmp;
  double **P;
  double *D;
  double **alpha;
  double *ptr1, *ptr2;
  double *laj_ptr;
  double dot_prod;
  double coef;
  double tmp_val1, tmp_val2;
  double maxabs;

  // Debug print
  if (debug_)
    ExEnv::out0() << indent << "Entered cphf" << endl;
  // End of debug print

  ////////////////////////////////////////////////////////////
  // Allocate and initialize various variables
  ////////////////////////////////////////////////////////////

  AP_matrix_tot = new double*[1];
  AP_matrix_tot[0] = new double[dimP];

  alpha = new double*[1];
  alpha[0] = new double[1];

  P = new double*[1];
  P[0] = new double[dimP];

  D = new double[nbasis*nbasis];

  // NB: Elements in Laj are ordered as (j*nvir + a)
  // since this ordering has been used with benefit in
  // MP2 gradient program
  ptr1 = P[0];
  ptr2 = P_sum_old;
  for (a=0; a<nvir; a++) {
    laj_ptr = &Laj[a];
    for (j=0; j<nocc; j++) {
      *ptr1++ = *laj_ptr/(eigval[a+nocc]-eigval[j]);
      *ptr2++ = 0.0;
      laj_ptr += nvir;
      }
    }

  for (i=0; i<nbasis; i++) {
    for (j=0; j<noso; j++) {
      if (j<nocc) Co(i,j) = scf_vector[i][j];
      else Cv(i,j-nocc) = scf_vector[i][j];
      }
    }

  /////////////////////////////////////////////////////////////////
  // Solve the CPHF equations (iteratively, with DIIS like method) 
  /////////////////////////////////////////////////////////////////

  i = 0;
  niter = 0;

  const int maxiter = 30;
  const int warniter = 1;
  while (niter < maxiter) { // Allow max maxiter iterations

    niter++;
    i++;
    if (debug_)
      ExEnv::out0() << indent << scprintf("niter: %i\n", niter);

    // First expand AP_matrix_tot, alpha and P with an extra row

    AP_matrix_tmp = new double *[i+1];
    if (!AP_matrix_tmp) {
      ExEnv::errn() << "Could not allocate AP_matrix_tmp" << endl;
      abort();
      }

    alpha_tmp = new double *[i+1];

    if (!alpha_tmp) {
      ExEnv::errn() << "Could not allocate alpha_tmp" << endl;
      abort();
      }

    P_tmp = new double *[i+1];

    if (!P_tmp) {
      ExEnv::errn() << "Could not allocate P_tmp" << endl;
      abort();
      }

    for (j=0; j<i; j++) {
      AP_matrix_tmp[j] = AP_matrix_tot[j];
      alpha_tmp[j] = alpha[j];
      P_tmp[j] = P[j];
      }

    AP_matrix_tmp[i] = new double[dimP];
    if (!AP_matrix_tmp[i]) {
      ExEnv::errn() << scprintf("Could not allocate AP_matrix_tmp[i], i = %i",i) << endl;
      abort();
      }
    delete[] AP_matrix_tot;
    AP_matrix_tot = AP_matrix_tmp;

    alpha_tmp[i] = new double[1];
    if (!alpha_tmp[i]) {
      ExEnv::errn() << scprintf("Could not allocate alpha_tmp[i], i = %i",i) << endl;
      abort();
      }
    delete[] alpha;
    alpha = alpha_tmp;

    P_tmp[i] = new double[dimP];
    if (!P_tmp[i]) {
      ExEnv::errn() << scprintf("Could not allocate P_tmp[i], i = %i",i) << endl;
      abort();
      }
    delete[] P;
    P = P_tmp;

    // Initialize P[i]
    for (j=0; j<dimP; j++) P[i][j] = 0.0;

    // Compute A*P[i-1] (called AP_matrix) which is required to get P[i]
    // A*P[i-1] is treated as a matrix to facilitate its computation
    // A*P[i-1] is put into row i-1 of AP_matrix_tot

    ptr1 = P[i-1];
    for (j=0; j<nvir; j++) {
      for (k=0; k<nocc; k++) {
        P_matrix->set_element(j,k,*ptr1++);  // Convert P[i-1] to RefSCMatrix
        }
      }
    D_matrix = Cv*P_matrix*Co.t();
#if 0
    D_matrix = D_matrix + D_matrix.t();
    D_matrix->convert(D);  // Convert D_matrix to double* D
    make_cs_gmat(G, D);
#else
    RefSymmSCMatrix sD(D_matrix.rowdim(), kit);
    sD.assign(0.0);
    sD.accumulate_symmetric_sum(D_matrix);
    make_cs_gmat_new(G, sD);
#endif
    AP_matrix = 2.0*Cv.t()*G*Co;

    ptr1 = AP_matrix_tot[i-1];
    for (j=0; j<nvir; j++) {
      for (k=0; k<nocc; k++) {
        tmp_val1 = AP_matrix->get_element(j,k)/(eigval[k]-eigval[j+nocc]);
        AP_matrix->set_element(j,k,tmp_val1);
        *ptr1++ = tmp_val1;
        }
      }
    // End of AP_matrix computation

    // Compute coefficients alpha[0],...,alpha[i-1]
    compute_alpha(i, AP_matrix_tot, alpha, P, eigval, nocc, nvir);

    // Compute the vector P_sum_new = alpha[0]P[0]+...+alpha[i-1]P[i-1]
    ptr1 = P_sum_new;
    for (j=0; j<dimP; j++) *ptr1++ = 0.0;
    for (j=0; j<i; j++) {
      tmp_val1 = alpha[j][0];
      ptr1 = P_sum_new;
      ptr2 = P[j];
      for (k=0; k<dimP; k++) {
        *ptr1++ += tmp_val1 * *ptr2++;
        }
      }


    /////////////////////////////////////////////////////////////
    // Test for convergence 
    // (based on RMS(P2aj_new - P2aj_old) 
    //  and max abs. val. of element)
    /////////////////////////////////////////////////////////////
    ptr1 = P_sum_new;
    ptr2 = P_sum_old;
    tmp_val1 = 0.0;
    maxabs = 0.0;
    for (j=0; j<dimP; j++) {
       tmp_val2 = *ptr1++ - *ptr2++;
       tmp_val1 += tmp_val2*tmp_val2;
       if (fabs(tmp_val2) > maxabs) maxabs = fabs(tmp_val2);
       }
    if (debug_) {
      ExEnv::out0() << indent << scprintf("RMS(P2aj_new-P2aj_old) = %12.10lf",
                                          sqrt((tmp_val1)/dimP))
           << endl;
      ExEnv::out0() << indent
           << scprintf("max. abs. element of (P2aj_new-P2aj_old) = %12.10lf",
                       maxabs)
           << endl;
      }
    if (sqrt(tmp_val1)/dimP < epsilon && maxabs < epsilon) break; // Converged

    // Put P_sum_new into P_sum old
    ptr1 = P_sum_new;
    ptr2 = P_sum_old;
    for (j=0; j<dimP; j++) {
      *ptr2++ = *ptr1++;
      }

    // Compute projection of A*P[i-1] on P[0],...,P[i-1]

    ptr1 = projctn;
    for (j=0; j<dimP; j++) *ptr1++ = 0.0;
    for (j=0; j<i; j++) {
      dot_prod = 0.0;
      ptr1 = P[j];
      for (k=0; k<dimP; k++) {
        tmp_val1 = *ptr1++;
        dot_prod += tmp_val1*tmp_val1; // Compute dot product <P[j]|P[j]>
        }
      ptr1 = P[j];
      coef = 0.0;
      for (k=0; k<nvir; k++) {
        for (l=0; l<nocc; l++) {
          coef += *ptr1++ * AP_matrix->get_element(k,l);
          }
        }
      coef /= dot_prod;
      ptr1 = P[j];
      ptr2 = projctn;
      for (k=0; k<dimP; k++) {
        *ptr2++ += coef * *ptr1++;
        }
      }

    // Compute P[i] as A*P[i-1] - projctn

    ptr1 = P[i];
    ptr2 = projctn;
    for (j=0; j<nvir; j++) {
      for (k=0; k<nocc; k++) {
        *ptr1++ = AP_matrix->get_element(j,k) - *ptr2++;
        }
      }

    /////////////////////////////////////////////
    // Test for convergence (based on norm(P[i])
    /////////////////////////////////////////////
    tmp_val1 = 0.0;
    for (l=0; l<dimP; l++) {
      tmp_val1 += P[niter][l]*P[niter][l];
      }
    tmp_val1 = sqrt(tmp_val1/dimP);
    if (debug_)
      ExEnv::out0() << indent
           << scprintf("norm(P[niter]) = %12.10lf", tmp_val1) << endl;
    if (tmp_val1 < epsilon) {  // Converged (if norm of new vector is zero)
      ExEnv::out0() << indent
                   << scprintf("CPHF: iter = %2d rms(P) = %12.10f eps = %12.10f",
                               niter, tmp_val1, epsilon)
                   << endl << endl;
      break;
      }

    if (niter >= warniter) {
      ExEnv::out0() << indent
           << scprintf("CPHF: iter = %2d rms(P) = %12.10f eps = %12.10f",
                       niter, tmp_val1, epsilon)
           << endl;
      }

    }

  ///////////////////////////////////////////////////////////////
  // If CPHF equations did not converge, exit with error message
  ///////////////////////////////////////////////////////////////
  if (niter == maxiter) {
    ExEnv::out0() << indent
         << "CPHF equations did not converge in " << maxiter << " iterations"
         << endl;
    abort();
    }

  /////////////////////////////////////////////////////
  // The converged vector is in P_sum_new;
  // Put elements into P2aj
  // NB: Elements in P2aj are ordered as (a*nocc + j);
  /////////////////////////////////////////////////////
  ptr1 = P_sum_new;
  for (i=0; i<nvir; i++) {
    for (j=0; j<nocc; j++) {
      P2aj->set_element(i,j,*ptr1++);
      }
    }

  // Debug print
  if (debug_)
    ExEnv::out0() << indent << "Exiting cphf" << endl;
  // End of debug print

  // Deallocate various arrays
  delete[] D;

  for (i=0; i<niter+1; i++) {
    delete[] AP_matrix_tot[i];
    delete[] alpha[i];
    delete[] P[i];
    }
  delete[] AP_matrix_tot;
  delete[] alpha;
  delete[] P;
  delete[] projctn;
  delete[] P_sum_new;
  delete[] P_sum_old;
}


static void
compute_alpha(int dim, double **AP, double **alpha,
              double **P, double *eigval, int nocc, int nvir)
{
  //////////////////////////////////////////////////////
  // Solve the linear system of equations C*X = B
  // where C is RefSCMatrix and X and B are RefSCVector
  // Put result (X) into array alpha
  //////////////////////////////////////////////////////
  int i, j, k;
  int vect_dim = nocc*nvir;

  double tmp1, tmp2;
  double *ptr1, *ptr2;
  double *norm = new double[dim]; // contains norms of vectors P[i], i=0,dim

  Ref<SCMatrixKit> kit = SCMatrixKit::default_matrixkit();
  RefSCDimension C_dim(new SCDimension(dim));

  RefSCMatrix C(C_dim,C_dim,kit);
  RefSCVector B(C_dim,kit);
  RefSCVector X(C_dim,kit);

  // Compute norms of vectors P[i] and put into norm
  for (i=0; i<dim; i++) norm[i] = 0.0;
  ptr1 = norm;
  for (i=0; i<dim; i++) {
    ptr2 = P[i];
    for (j=0; j<vect_dim; j++) {
      *ptr1 += *ptr2 * *ptr2;
       ptr2++;
       }
    ptr1++;
    }
  for (i=0; i<dim; i++) norm[i] = sqrt(norm[i]);
    
  // Construct matrix C
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      tmp1 = 0.0;
      ptr1 = P[i];
      ptr2 = AP[j];
      for (k=0; k<vect_dim; k++) {
        tmp1 -= *ptr1++ * *ptr2++;
        }
      if (i == j) {
        ptr1 = P[i];
        for (k=0; k<vect_dim; k++) {
          tmp2 = *ptr1++;
          tmp1 += tmp2*tmp2;
          }
        }
      C->set_element(i,j,tmp1/(norm[i]*norm[j]));
      }
    }

  // Construct vector B
  B->set_element(0,norm[0]);
  for (i=1; i<dim; i++) B->set_element(i,0.0);

  // Compute X = inv(C)*B
  X = C.i()*B;

  // Put elements of X into alpha
  for (i=0; i<dim; i++) {
    alpha[i][0] = X->get_element(i)/norm[i];
    }

  delete[] norm;
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
