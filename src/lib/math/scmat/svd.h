//
// svd.h
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

#ifndef _chemistry_qc_mbptr12_svd_h
#define _chemistry_qc_mbptr12_svd_h

#include <math/scmat/matrix.h>
#include <math/scmat/blas.h>

namespace sc {

  /** Uses LAPACK's DGESVD to perform SVD of A:
   A = U * Sigma * V
   */
  void lapack_svd(const RefSCMatrix& A, RefSCMatrix& U, RefDiagSCMatrix& Sigma,
                  RefSCMatrix& V);

  /** Uses LAPACK's DSPSVX to solve symmetric non-definite linear system AX = B, where
      B is a single vector
      @return The estimate of the reciprocal condition number of the matrix A.
   */
  double lapack_linsolv_symmnondef(const RefSymmSCMatrix& A, RefSCVector& X,
                                   const RefSCVector& B);

  /** Uses LAPACK's DSPSVX to solve symmetric non-definite linear system AX = B, where
   B is a set of vectors
      @return The estimate of the reciprocal condition number of the matrix A.
   */
  double lapack_linsolv_symmnondef(const RefSymmSCMatrix& A, RefSCMatrix& X,
                                   const RefSCMatrix& B);

  /** Same as above, except uses C-style arrays already.
   A is in packed upper-triangular form, nA is the dimension of A,
   X and B are matrices with ncolB rows and nA columns.
      @return The estimate of the reciprocal condition number of the matrix A.
   */
  double lapack_linsolv_symmnondef(const double* AP, int nA, double* Xt,
                                   const double* Bt, int ncolB);

  /** invert symmetric non-definite matrix using DSPTRF LAPACK routine
   that implements the Bunch-Kaufman diagonal pivoting method.

   \param condition_number_threshold when the estimate of the condition number (see lapack function DSPCON)
   exceeds this threshold print a warning to ExEnv::err0(). No warning will be produced if the threshold is 0.0.
   Negative threshold will prompt the estimate of the condition number of A to be printed out to ExEnv::out0().
   */
  void lapack_invert_symmnondef(RefSymmSCMatrix& A,
                                double condition_number_threshold = 0.0);

  /** Compute factorization of a symmetric non-definite matrix using DSPTRF LAPACK routine
   that implements the Bunch-Kaufman diagonal pivoting method.
   The resulting factorization (AF and ipiv) can be used to solve a linear system A X = B
   using lapack_linsolv_dpf_symmnondef();

   \param condition_number_threshold when the estimate of the condition number (see lapack function DSPCON)
   exceeds this threshold print a warning to ExEnv::err0(). No warning will be produced if the threshold is 0.0.
   Negative threshold will prompt the estimate of the condition number of A to be printed out to ExEnv::out0().

   \param AF array of size A.dim().n() * (A.dim().n() + 1)/2. On output contains the desired factorization.

   \param ipiv integer array of sizeA.dim().n()

   */
  void lapack_dpf_symmnondef(const RefSymmSCMatrix& A, double* AF, blasint* ipiv,
                             double condition_number_threshold = 0.0);

  /** invert symmetric positive-definite matrix using DPPTRF LAPACK routine
   that implements the Cholesky method.

   \param condition_number_threshold when the estimate of the condition number (see lapack function DPPCON)
   exceeds this threshold print a warning to ExEnv::err0(). No warning will be produced if the threshold is 0.0.
   Negative threshold will prompt the estimate of the condition number of A to be printed out to ExEnv::out0().
   */
  void lapack_invert_symmposdef(RefSymmSCMatrix& A,
                                double condition_number_threshold = 0.0);

  /**
   * Solves a symmetric indefinite system of linear equations AX=B,
   * where A is held in packed storage, using the factorization computed by lapack_dpf_symmnondef
   * @param A
   * @param nA
   * @param AF
   * @param ipiv
   * @param Xt
   * @param Bt
   * @param ncolB
   * @param refine set to false to avoid the iterative refinement of the solution that was obtained by substitution (using LAPACK dsptrs function)
   */
  void lapack_linsolv_dpf_symmnondef(const double* A, int nA, const double* AF,
                                     const blasint* ipiv, double* Xt,
                                     const double* Bt, int ncolB,
                                     bool refine = true);

  /** Solves symmetric non-definite linear system AX = B, where B is a RefSCVector, using Jacobi solver.
   *  Note that Jacobi solver rarely converges, hence it should not be used.
   *
   *       @return The estimate of the reciprocal condition number of the matrix A.
   */
  double linsolv_symmnondef_jacobi(const RefSymmSCMatrix& A, RefSCVector& X,
                                   const RefSCVector& B);

  /** Solves symmetric non-definite linear system AX = B, where B is a RefSCVector, using conjugate gradient solver
   *
   *       @return The estimate of the reciprocal condition number of the matrix A.
   */
  double linsolv_symmnondef_cg(const RefSymmSCMatrix& A, RefSCVector& X,
                               const RefSCVector& B);

  /** Compute factorization of a symmetric positive-definite matrix using DPPTRF LAPACK routine
   that implements the Cholesky method.
   The resulting factorization (AF and ipiv) can be used to solve a linear system A X = B
   using lapack_linsolv_cholesky_symmposdef();

   \param condition_number_threshold when the estimate of the condition number (see lapack function DSPCON)
   exceeds this threshold print a warning to ExEnv::err0(). No warning will be produced if the threshold is 0.0.
   Negative threshold will prompt the estimate of the condition number of A to be printed out to ExEnv::out0().

   \param AF array of size A.dim().n() * (A.dim().n() + 1)/2. On output contains the desired factorization.

   */
  void lapack_cholesky_symmposdef(const RefSymmSCMatrix& A, double* AF,
                                  double condition_number_threshold = 0.0);

  /**
   * Solves a symmetric indefinite system of linear equations AX=B,
   * where A is held in packed storage, using the factorization computed by lapack_cholesky_symmnondef
   * @param A
   * @param nA
   * @param AF
   * @param Xt
   * @param Bt
   * @param ncolB
   * @param refine set to false to avoid the iterative refinement of the solution that was obtained by substitution (using LAPACK dsptrs function)
   */
  void lapack_linsolv_cholesky_symmposdef(const double* A, int nA,
                                          const double* AF, double* Xt,
                                          const double* Bt, int ncolB,
                                          bool refine = true);

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:


