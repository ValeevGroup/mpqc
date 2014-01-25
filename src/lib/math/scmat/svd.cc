//
// svd.cc
//
// Copyright (C) 2005 Edward Valeev
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

#include <cmath>
#include <cassert>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <util/misc/formio.h>
#include <util/misc/consumableresources.h>
#include <util/misc/scexception.h>
#include <math/scmat/blas.h>
#include <math/scmat/lapack.h>
#include <math/scmat/svd.h>

using namespace std;
using namespace sc;

namespace sc {

  /** Uses LAPACK's DGESVD to perform SVD of A:
   A = U * Sigma * V
   */
  void lapack_svd(const RefSCMatrix& A, RefSCMatrix& U, RefDiagSCMatrix& Sigma,
                  RefSCMatrix& V) {
    /*
     In terms of C matrixes, Fortran77 call is to compute SVD of A^t

     A^t = V^t * Sigma^t * U^t

     DGESVD uses notation A = U * Sigma * Vt, hence the C-F77 correspondence is as follows:

     C     F77
     ---    ---
     A      A^t
     V      U
     Sigma  Sigma
     U      Vt

     */
    blasint m = A.nrow();
    blasint n = A.ncol();
    blasint nsigma = m < n ? m : n;
    blasint lwork = m < n ? 5 * n : 5 * m;
    char jobvt = U.nonnull() ? 'A' : 'N';
    char jobu = V.nonnull() ? 'A' : 'N';

    double* a_data = allocate<double>(m * n);
    A.convert(a_data);

    double* u_data = 0;
    if (U.nonnull()) {
      u_data = allocate<double>(m * m);
    }
    double* v_data = 0;
    if (V.nonnull()) {
      v_data = allocate<double>(n * n);
    }
    double* s_data = new double[nsigma];
    double *work = new double[lwork];
    blasint info = 0;

    F77_DGESVD(&jobu, &jobvt, &n, &m, a_data, &n, s_data, v_data, &n, u_data,
               &m, work, &lwork, &info);

    if (info != 0) {
      delete[] work;
      deallocate(a_data);
      delete[] s_data;
      if (u_data)
        deallocate(u_data);
      if (v_data)
        deallocate(v_data);
      throw std::runtime_error("lapack_svd -- call to LAPACK routine failed");
    }

    Sigma.assign(s_data);
    if (U.nonnull()) {
      U.assign(u_data);
    }
    if (V.nonnull()) {
      V.assign(v_data);
    }

    delete[] work;
    deallocate(a_data);
    delete[] s_data;
    if (u_data)
      deallocate(u_data);
    if (v_data)
      deallocate(v_data);
  }

  /** Uses LAPACK's DSPSVX to solve symmetric non-definite linear system AX = B, where B is a RefSCVector
   */
  double lapack_linsolv_symmnondef(const RefSymmSCMatrix& A, RefSCVector& X,
                                   const RefSCVector& B) {
    blasint n = A.n();
    if (n != B.n())
      throw std::runtime_error(
                               "lapack_linsolv_symmnondef() -- dimensions of A and B do not match");
    if (n != X.n())
      throw std::runtime_error(
                               "lapack_linsolv_symmnondef() -- dimensions of A and X do not match");

    // convert A to packed upper triangular form
    blasint ntri = n * (n + 1) / 2;
    double* AP = allocate<double>(ntri);
    A.convert(AP);
    double* AFP = allocate<double>(ntri);
    blasint* ipiv = new blasint[n];
    double* BB = new double[n];
    B.convert(BB);
    double* XX = new double[n];

    char fact = 'N';
    char uplo = 'U';
    blasint nrhs = 1;
    double rcond = 0.0;
    double* ferr = new double[1];
    double* berr = new double[1];
    double* work = new double[3 * n];
    blasint* iwork = new blasint[n];
    blasint info = 0;

    F77_DSPSVX(&fact, &uplo, &n, &nrhs, AP, AFP, ipiv, BB, &n, XX, &n, &rcond,
               ferr, berr, work, iwork, &info);

    if (info) {

      delete[] ferr;
      delete[] berr;
      delete[] work;
      delete[] iwork;

      deallocate(AP);
      deallocate(AFP);
      delete[] BB;
      delete[] XX;

      if (info > 0 && info <= n)
        throw std::runtime_error(
                                 "lapack_linsolv_symmnondef() -- matrix A has factors which are exactly zero");
      if (info == n + 1)
        throw std::runtime_error(
                                 "lapack_linsolv_symmnondef() -- matrix A is singular");
      if (info < 0)
        throw std::runtime_error(
                                 "lapack_linsolv_symmnondef() -- one of the arguments to F77_DSPSVX is invalid");
    } else {
      X.assign(XX);

      delete[] ferr;
      delete[] berr;
      delete[] work;
      delete[] iwork;

      deallocate(AP);
      deallocate(AFP);
      delete[] BB;
      delete[] XX;
    }
    return rcond;
  }

  /** Uses LAPACK's DSPSVX to solve symmetric non-definite linear system AX = B, where B is a RefSCMatrix.
   Dimensions of X and B must match.
   */
  double lapack_linsolv_symmnondef(const RefSymmSCMatrix& A, RefSCMatrix& X,
                                   const RefSCMatrix& B) {
    blasint n = A.n();
    blasint nrhs = B.ncol();
    if (n != B.nrow())
      throw std::runtime_error(
                               "lapack_linsolv_symmnondef() -- dimensions of A and B do not match");
    if (n != X.nrow())
      throw std::runtime_error(
                               "lapack_linsolv_symmnondef() -- dimensions of A and X do not match");
    if (X.ncol() != nrhs)
      throw std::runtime_error(
                               "lapack_linsolv_symmnondef() -- dimensions of B and X do not match");

    // convert A to packed upper triangular form
    blasint ntri = n * (n + 1) / 2;
    double* AP = allocate<double>(ntri);
    A.convert(AP);
    double* AFP = allocate<double>(ntri);
    blasint* ipiv = new blasint[n];
    double* BB = allocate<double>(nrhs * n);
    B.t().convert(BB);
    double* XX = allocate<double>(nrhs * n);

    char fact = 'N';
    char uplo = 'U';
    double rcond = 0.0;
    double* ferr = new double[nrhs];
    double* berr = new double[nrhs];
    double* work = new double[3 * n];
    blasint* iwork = new blasint[n];
    blasint info = 0;

    F77_DSPSVX(&fact, &uplo, &n, &nrhs, AP, AFP, ipiv, BB, &n, XX, &n, &rcond,
               ferr, berr, work, iwork, &info);

    if (info) {

      delete[] ferr;
      delete[] berr;
      delete[] work;
      delete[] iwork;

      deallocate(AP);
      deallocate(AFP);
      deallocate(BB);
      deallocate(XX);

      if (info > 0 && info <= n)
        throw std::runtime_error(
                                 "lapack_linsolv_symmnondef() -- matrix A has factors which are exactly zero");
      if (info == n + 1)
        throw std::runtime_error(
                                 "lapack_linsolv_symmnondef() -- matrix A is singular");
      if (info < 0)
        throw std::runtime_error(
                                 "lapack_linsolv_symmnondef() -- one of the arguments to F77_DSPSVX is invalid");
    } else {
      RefSCMatrix Xt = X.t();
      Xt.assign(XX);
      X = Xt.t();
      Xt = 0;

      delete[] ferr;
      delete[] berr;
      delete[] work;
      delete[] iwork;

      deallocate(AP);
      deallocate(AFP);
      deallocate(BB);
      deallocate(XX);
    }
    return rcond;
  }

  /** Same as above, except uses C-style arrays already.
   A is in packed upper-triangular form, nA is the dimension of A,
   X and B are matrices with ncolB rows and nA columns.
   */
  double lapack_linsolv_symmnondef(const double* AP, int nA, double* Xt,
                                   const double* Bt, int ncolB) {
    const blasint n = nA;
    const blasint nrhs = ncolB;

    const blasint ntri = n * (n + 1) / 2;
    // AFP is a scratch array
    double* AFP = allocate<double>(ntri);
    blasint* ipiv = new blasint[n];

    char fact = 'N';
    char uplo = 'U';
    double rcond = 0.0;
    double* ferr = new double[nrhs];
    double* berr = new double[nrhs];
    double* work = new double[3 * n];
    blasint* iwork = new blasint[n];
    blasint info = 0;

    F77_DSPSVX(&fact, &uplo, &n, &nrhs, AP, AFP, ipiv, Bt, &n, Xt, &n, &rcond,
               ferr, berr, work, iwork, &info);

    delete[] ferr;
    delete[] berr;
    delete[] work;
    delete[] iwork;

    deallocate(AFP);

    if (info) {

      if (info > 0 && info <= n)
        throw std::runtime_error(
                                 "lapack_linsolv_symmnondef() -- matrix A has factors which are exactly zero");
      if (info == n + 1)
        throw std::runtime_error(
                                 "lapack_linsolv_symmnondef() -- matrix A is singular");
      if (info < 0)
        throw std::runtime_error(
                                 "lapack_linsolv_symmnondef() -- one of the arguments to F77_DSPSVX is invalid");
    }

    return rcond;
  }

  /** Solves symmetric non-definite linear system AX = B, where B is a RefSCVector, using Jacobi solver
   *  This is likely to diverge
   */
  double linsolv_symmnondef_jacobi(const RefSymmSCMatrix& A, RefSCVector& X,
                                   const RefSCVector& B) {
    const blasint n = A.n();
    if (n != B.n())
      throw std::runtime_error(
                               "linsolv_symmnondef_jacobi() -- dimensions of A and B do not match");
    if (n != X.n())
      throw std::runtime_error(
                               "linsolv_symmnondef_jacobi() -- dimensions of A and X do not match");

    // convert A to packed lower triangular form (upper triangular in fortran)
    const blasint ntri = n * (n + 1) / 2;
    double* RP = allocate<double>(ntri);
    A.convert(RP);
    for(blasint i=0, ii=0; i<n; ii+=(++i + 1))
      RP[ii] = 0.0;
    // convert B to a dense vector
    double* BB = allocate<double>(n);
    B.convert(BB);

    // soluton vectors
    double* XX_i = allocate<double>(n);
    double* XX_ip1 = allocate<double>(n);

    // extract the diagonal inverse
    double* DDinv = allocate<double>(n);
    for(blasint i=0; i<n; ++i)
      DDinv[i] = 1.0 / A(i,i);

    // approximate the condition number as the ratio of the min and max elements of the diagonal
    const double DDinv_min = *std::min_element(DDinv, DDinv+n);
    const double DDinv_max = *std::max_element(DDinv, DDinv+n);
    const double DD_cond = DDinv_max / DDinv_min;
    const double conv_threshold = 1e-15 * DD_cond; // estimate of how tightly the system can be converged

    char uplo = 'U';
    const double minusone = -1.0;
    const double one = 1.0;
    const blasint ione = 1;
    bool converged = false;
    const unsigned int max_niter = 50000;

    // starting guess: x_0 = D^-1 . b
    for(blasint i=0; i<n; ++i)
      XX_i[i] = BB[i] * DDinv[i];
    // exact solution
//    lapack_linsolv_symmnondef(A, X, B);
//    for(int i=0; i<n; ++i)
//      XX_i[i] = X(i);

    unsigned int iter = 0;
    double norm = 0.0;
    while (not converged) {

      // XX_ip1 = b - R . x_i
      std::copy(BB, BB+n, XX_ip1);
      F77_DSPMV(&uplo, &n, &minusone, RP, XX_i, &ione, &one, XX_ip1, &ione);

      // x_i+1 = D^-1 . XX_ip1 = D^-1 . (b - R . x_i)
      for(blasint i=0; i<n; ++i)
        XX_ip1[i] *= DDinv[i];

      norm = 0.0;
      for(blasint i=0; i<n; ++i) {
        const double dx = (XX_ip1[i] - XX_i[i]);
        norm += dx*dx;
      }
      norm = sqrt(norm) / n; // divide by n to make it dimension-independent
      if (norm < conv_threshold)
        converged = true;
      ++iter;
      //ExEnv::out0() << indent << scprintf("iter=%d dnorm=%25.15lf",iter,norm) << std::endl;
      if (iter >= max_niter) {
        deallocate(RP);
        deallocate(BB);
        deallocate(DDinv);
        deallocate(XX_i);
        deallocate(XX_ip1);
        throw MaxIterExceeded("linsolv_symmnondef_jacobi() -- did not converge within the max # of iterations",
                              __FILE__, __LINE__,
                              max_niter);
      }

      std::swap(XX_ip1, XX_i);
    }

    X.assign(XX_ip1);

    deallocate(RP);
    deallocate(BB);
    deallocate(DDinv);
    deallocate(XX_i);
    deallocate(XX_ip1);

    return 1.0 / DD_cond;
  }

  /** Solves symmetric non-definite linear system AX = B, where B is a RefSCVector, using conjugate gradient solver
   */
  double linsolv_symmnondef_cg(const RefSymmSCMatrix& A, RefSCVector& X,
                                   const RefSCVector& B) {
    const blasint n = A.n();
    if (n != B.n())
      throw std::runtime_error(
                               "linsolv_symmnondef_cg() -- dimensions of A and B do not match");
    if (n != X.n())
      throw std::runtime_error(
                               "linsolv_symmnondef_cg() -- dimensions of A and X do not match");

    // convert A to packed lower triangular form (upper triangular in fortran)
    const blasint ntri = n * (n + 1) / 2;
    double* AP = allocate<double>(ntri);
    A.convert(AP);
    // convert B to a dense vector
    double* BB = allocate<double>(n);
    B.convert(BB);

    // soluton vector
    double* XX_i = allocate<double>(n);
    // residual vector
    double* RR_i = allocate<double>(n);
    // preconditioned residual vector
    double* ZZ_i = allocate<double>(n);
    // direction vector
    double* PP_i = allocate<double>(n);
    double* APP_i = allocate<double>(n);

    // extract the diagonal |inverse|
    double* DDinv = allocate<double>(n);
    for(blasint i=0; i<n; ++i)
      DDinv[i] = 1.0 / fabs(A(i,i));
    // approximate the condition number as the ratio of the min and max elements of the diagonal
    const double DDinv_min = *std::min_element(DDinv, DDinv+n);
    const double DDinv_max = *std::max_element(DDinv, DDinv+n);
    const double DD_cond = DDinv_max / DDinv_min;
    const double conv_threshold = 1e-15 * DD_cond; // estimate of how tightly the system can be converged

    char uplo = 'U';
    const double minusone = -1.0;
    const double one = 1.0;
    const double zero = 0.0;
    const blasint ione = 1;
    bool converged = false;
    const unsigned int max_niter = 50000;

    // starting guess: x_0 = D^-1 . b
    std::transform(BB, BB+n, DDinv, XX_i, std::multiplies<double>());
    // exact solution
//    lapack_linsolv_symmnondef(A, X, B);
//    for(int i=0; i<n; ++i)
//      XX_i[i] = X(i);

    // r_0 = b - A . x_0
    std::copy(BB, BB+n, RR_i);
    F77_DSPMV(&uplo, &n, &minusone, AP, XX_i, &ione, &one, RR_i, &ione);

    // z_0 = D^-1 . r_0
    std::transform(RR_i, RR_i+n, DDinv, ZZ_i, std::multiplies<double>());

    // p_0 = z_0
    std::copy(ZZ_i, ZZ_i+n, PP_i);

    unsigned int iter = 0;
    while (not converged) {

      // alpha_i = (r_i . z_i) / (p_i . A . p_i)
      const double rz_norm2 = F77_DDOT(&n, RR_i, &ione,
                                       ZZ_i, &ione);
      F77_DSPMV(&uplo, &n, &one, AP, PP_i, &ione, &zero, APP_i, &ione);
      const double pAp_i = F77_DDOT(&n, PP_i, &ione,
                                    APP_i, &ione);
      const double alpha_i = rz_norm2 / pAp_i;

      // x_i += alpha_i p_i
      F77_DAXPY(&n, &alpha_i, PP_i, &ione, XX_i, &ione);

      // r_i -= alpha_i Ap_i
      const double minus_alpha_i = -alpha_i;
      F77_DAXPY(&n, &minus_alpha_i, APP_i, &ione, RR_i, &ione);

      const double r_ip1_norm = sqrt(F77_DDOT(&n, RR_i, &ione,
                                              RR_i, &ione)) / n;  // divide by n to make this dimension-independent
      if (r_ip1_norm < conv_threshold)
        converged = true;

      // z_i = D^-1 . r_i
      std::transform(RR_i, RR_i+n, DDinv, ZZ_i, std::multiplies<double>());
      const double rz_ip1_norm2 = F77_DDOT(&n, RR_i, &ione,
                                               ZZ_i, &ione);

      const double beta_i = rz_ip1_norm2 / rz_norm2;

      // p_i = z_i+1 + beta_i p_i
      // 1) scale by beta_i
      // 2) add z_i+1 (i.e. current contents of z_i)
      F77_DSCAL(&n, &beta_i, PP_i, &ione);
      F77_DAXPY(&n, &one, ZZ_i, &ione, PP_i, &ione);

      ++iter;
      ExEnv::out0() << indent << scprintf("iter=%d dnorm=%25.15lf",iter,r_ip1_norm) << std::endl;

      if (iter >= max_niter) {
        deallocate(AP);
        deallocate(BB);
        deallocate(XX_i);
        deallocate(RR_i);
        deallocate(PP_i);
        deallocate(ZZ_i);
        deallocate(APP_i);
        deallocate(DDinv);
        throw MaxIterExceeded("linsolv_symmnondef_cg() -- did not converge within the max # of iterations",
                              __FILE__, __LINE__,
                              max_niter);
      }
    } // solver loop

    X.assign(XX_i);

    deallocate(AP);
    deallocate(BB);
    deallocate(XX_i);
    deallocate(RR_i);
    deallocate(PP_i);
    deallocate(ZZ_i);
    deallocate(APP_i);
    deallocate(DDinv);

    return 1.0 / DD_cond;
  }

  /** Computed eigenvalues of A and returns how many are below threshold. Uses LAPACK's DSYEVD.
   */
  int lapack_evals_less_eps(const RefSymmSCMatrix& A, double threshold) {
    blasint n = A.n();
    // convert A to square form
    double* Asq = allocate<double>(n * n);
    double** Arows = new double*[n];
    blasint ioff = 0;
    for (blasint i = 0; i < n; i++, ioff += n)
      Arows[i] = Asq + ioff;
    A.convert(Arows);
    delete[] Arows;
    double* evals = new double[n];

    char need_evecs = 'N';
    char uplo = 'U';
    const blasint lwork = 2 * n + 1;
    double* work = new double[lwork];
    const blasint liwork = 1;
    blasint* iwork = new blasint[liwork];
    blasint info = 0;

    F77_DSYEVD(&need_evecs, &uplo, &n, Asq, &n, evals, work, &lwork, iwork,
               &liwork, &info);

    if (info) {

      delete[] work;
      delete[] iwork;
      deallocate(Asq);
      delete[] evals;

      throw sc::ProgrammingError(
                                 "lapack_evals_less_eps() -- call to LAPACK's DSYEVD failed",
                                 __FILE__, __LINE__);
    }

    blasint nevals_lessthan_threshold = 0;
    for (blasint i = 0; i < n; i++) {
      if (evals[i] < threshold)
        nevals_lessthan_threshold++;
    }

    delete[] work;
    delete[] iwork;
    deallocate(Asq);
    delete[] evals;

    return nevals_lessthan_threshold;
  }

  void lapack_invert_symmnondef(RefSymmSCMatrix& A,
                                double condition_number_threshold) {
    const blasint n = A.dim().n();
    const blasint ntri = n * (n + 1) / 2;
    char uplo = 'U';
    std::vector<double> AF(ntri);
    std::vector<blasint> ipiv(n);
    blasint info;
    std::vector<double> work(2 * n);

    // compute factorization of A
    lapack_dpf_symmnondef(A, &(AF[0]), &(ipiv[0]), condition_number_threshold);

    // compute the inverse
    F77_DSPTRI(&uplo, &n, &(AF[0]), &(ipiv[0]), &(work[0]), &info);
    if (info) {
      if (info < 0)
        throw std::runtime_error(
                                 "lapack_invert_symmnondef() -- one of the arguments to F77_DSPTRI is invalid");
      if (info > 0)
        throw std::runtime_error(
                                 "lapack_invert_symmnondef() -- matrix A has factors which are exactly zero");
      MPQC_ASSERT(false); // unreachable
    }

    A.assign(&(AF[0]));
  }

  void lapack_dpf_symmnondef(const RefSymmSCMatrix& A, double* AF, blasint* ipiv,
                             double condition_number_threshold) {

    const blasint n = A.dim().n();
    char uplo = 'U';
    blasint info;
    std::vector<double> work(2 * n);
    std::vector<blasint> iwork(n);

    // compute the infinity-norm of A
    A.convert(AF);
    char norm = 'I';
    const double anorm = F77_DLANSP(&norm, &uplo, &n, AF, &(work[0]));

    // factorize A = U.D.Ut
    F77_DSPTRF(&uplo, &n, AF, ipiv, &info);
    if (info) {
      if (info < 0)
        throw std::runtime_error(
                                 "lapack_dpf_symmnondef() -- one of the arguments to F77_DSPTRF is invalid");
      if (info > 0)
        throw std::runtime_error(
                                 "lapack_dpf_symmnondef() -- matrix A has factors which are exactly zero");
      MPQC_ASSERT(false); // unreachable
    }

    // estimate the condition number
    double rcond;
    F77_DSPCON(&uplo, &n, AF, ipiv, &anorm, &rcond, &(work[0]), &(iwork[0]),
               &info);
    if (info) {
      if (info < 0)
        throw std::runtime_error(
                                 "lapack_dpf_symmnondef() -- one of the arguments to F77_DSPCON is invalid");
      MPQC_ASSERT(false); // unreachable
    }
    // if the condition number is above the threshold or its inverse below the working precision, throw
    if (condition_number_threshold != 0.0) {
      if (condition_number_threshold < 0.0)
        ExEnv::out0() << indent
            << "condition number estimate in lapack_dpf_symmnondef() = " << 1.0
            / rcond << std::endl;
      else if (1.0 / rcond > condition_number_threshold)
        ExEnv::err0() << indent
            << "WARNING: large condition number in lapack_dpf_symmnondef(): threshold = "
            << condition_number_threshold << " actual = " << 1.0 / rcond
            << std::endl;
      const char epsilon = 'E';
      if (rcond < F77_DLAMCH(&epsilon))
        ExEnv::err0() << indent
            << "WARNING: condition number in lapack_dpf_symmnondef() exceeds the working precision"
            << std::endl;
    }

  }

  void lapack_linsolv_dpf_symmnondef(const double* A, int nA, const double* AF,
                                     const blasint* ipiv, double* Xt,
                                     const double* Bt, int ncolB, bool refine) {
    const blasint n = nA;
    const char uplo = 'U';
    blasint info;
    const blasint ncolBf = ncolB;

    // copy B into X
    const char full = 'F';
    F77_DLACPY(&full, &n, &ncolBf, Bt, &n, Xt, &n, &info);

    // solve the linear system
    F77_DSPTRS(&uplo, &n, &ncolBf, AF, ipiv, Xt, &n, &info);

    if (refine) {
      // Use iterative refinement to improve the computed solutions and
      // compute error bounds and backward error estimates for them.
      std::vector<double> ferr(ncolB);
      std::vector<double> berr(ncolB);
      std::vector<double> work(3 * n);
      std::vector<blasint> iwork(n);
      F77_DSPRFS(&uplo, &n, &ncolBf, A, AF, ipiv, Bt, &n, Xt, &n, &(ferr[0]),
                 &(berr[0]), &(work[0]), &(iwork[0]), &info);
    }

  }

  void lapack_cholesky_symmposdef(const RefSymmSCMatrix& A, double* AF,
                                  double condition_number_threshold) {

    const blasint n = A.dim().n();
    char uplo = 'U';
    blasint info;
    std::vector<double> work(3 * n);
    std::vector<blasint> iwork(n);

    // compute the infinity-norm of A
    A.convert(AF);
    char norm = 'I';
    const double anorm = F77_DLANSP(&norm, &uplo, &n, AF, &(work[0]));

    // factorize A = Ut.U
    F77_DPPTRF(&uplo, &n, AF, &info);
    if (info) {
      if (info < 0)
        throw std::runtime_error(
                                 "lapack_cholesky_symmposdef() -- one of the arguments to F77_DPPTRF is invalid");
      if (info > 0)
        throw std::runtime_error(
                                 "lapack_cholesky_symmposdef() -- matrix A has factors which are negative");
      MPQC_ASSERT(false); // unreachable
    }

    // estimate the condition number
    double rcond;
    F77_DPPCON(&uplo, &n, AF, &anorm, &rcond, &(work[0]), &(iwork[0]), &info);
    if (info) {
      if (info < 0)
        throw std::runtime_error(
                                 "lapack_cholesky_symmposdef() -- one of the arguments to F77_DPPCON is invalid");
      MPQC_ASSERT(false); // unreachable
    }
    // if the condition number is above the threshold or its inverse below the working precision, throw
    if (condition_number_threshold != 0.0) {
      if (condition_number_threshold < 0.0)
        ExEnv::out0() << indent
            << "condition number estimate in lapack_cholesky_symmposdef() = "
            << 1.0 / rcond << std::endl;
      else if (1.0 / rcond > condition_number_threshold)
        ExEnv::err0() << indent
            << "WARNING: large condition number in lapack_cholesky_symmposdef(): threshold = "
            << condition_number_threshold << " actual = " << 1.0 / rcond
            << std::endl;
      const char epsilon = 'E';
      if (rcond < F77_DLAMCH(&epsilon))
        ExEnv::err0() << indent
            << "WARNING: condition number in lapack_cholesky_symmposdef() exceeds the working precision"
            << std::endl;
    }

  }

  void lapack_linsolv_cholesky_symmposdef(const double* A, int nA,
                                          const double* AF, double* Xt,
                                          const double* Bt, int ncolB,
                                          bool refine) {
    const blasint n = nA;
    const char uplo = 'U';
    blasint info;
    const blasint ncolBf = ncolB;

    // copy B into X
    const char full = 'F';
    F77_DLACPY(&full, &n, &ncolBf, Bt, &n, Xt, &n, &info);

    // solve the linear system
    F77_DPPTRS(&uplo, &n, &ncolBf, AF, Xt, &n, &info);

    if (refine) {
      // Use iterative refinement to improve the computed solutions and
      // compute error bounds and backward error estimates for them.
      std::vector<double> ferr(ncolB);
      std::vector<double> berr(ncolB);
      std::vector<double> work(3 * n);
      std::vector<blasint> iwork(n);
      F77_DPPRFS(&uplo, &n, &ncolBf, A, AF, Bt, &n, Xt, &n, &(ferr[0]),
                 &(berr[0]), &(work[0]), &(iwork[0]), &info);
    }

  }

  void lapack_invert_symmposdef(RefSymmSCMatrix& A,
                                double condition_number_threshold) {
    const blasint n = A.dim().n();
    const blasint ntri = n * (n + 1) / 2;
    char uplo = 'U';
    std::vector<double> AF(ntri);
    blasint info;

    // compute Cholesky factorization of A
    lapack_cholesky_symmposdef(A, &(AF[0]), condition_number_threshold);

    // compute the inverse
    F77_DPPTRI(&uplo, &n, &(AF[0]), &info);
    if (info) {
      if (info < 0)
        throw std::runtime_error(
                                 "lapack_invert_symmposdef() -- one of the arguments to F77_DPPTRI is invalid");
      if (info > 0)
        throw std::runtime_error(
                                 "lapack_invert_symmposdef() -- matrix A has factors which are exactly zero");
      MPQC_ASSERT(false); // unreachable
    }

    A.assign(&(AF[0]));
  }
};

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
