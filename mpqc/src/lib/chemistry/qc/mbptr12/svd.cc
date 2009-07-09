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

#ifdef __GNUG__
#pragma implementation
#endif

#include <cmath>
#include <cassert>
#include <stdexcept>
#include <util/misc/formio.h>
#include <util/class/scexception.h>
#include <chemistry/qc/mbptr12/lapack.h>
#include <chemistry/qc/mbptr12/svd.h>

using namespace std;
using namespace sc;

namespace sc {

  namespace exp {

    /** Uses LAPACK's DGESVD to perform SVD of A:
        A = U * Sigma * V
    */
    void lapack_svd(const RefSCMatrix& A,
                    RefSCMatrix& U,
                    RefDiagSCMatrix& Sigma,
                    RefSCMatrix& V)
    {
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
      int m = A.nrow();
      int n = A.ncol();
      int nsigma = m<n ? m : n;
      int lwork = m<n ? 5*n : 5*m;
      char jobvt = U ? 'A' : 'N';
      char jobu = V ? 'A' : 'N';

      double* a_data = new double[m*n];
      A.convert(a_data);

      double* u_data = 0;
      if (U) {
        u_data = new double[m*m];
      }
      double* v_data = 0;
      if (V) {
        v_data = new double[n*n];
      }
      double* s_data = new double[nsigma];
      double *work = new double[lwork];
      int info = 0;

      F77_DGESVD(&jobu, &jobvt, &n, &m,
                 a_data, &n, s_data,
                 v_data, &n,
                 u_data, &m,
                 work, &lwork, &info);

      if (info !=0) {
        delete[] work;
        delete[] a_data;
        delete[] s_data;
        if (u_data) delete[] u_data;
        if (v_data) delete[] v_data;
        throw std::runtime_error("lapack_svd -- call to LAPACK routine failed");
      }

      Sigma.assign(s_data);
      if (U) {
        U.assign(u_data);
      }
      if (V) {
        V.assign(v_data);
      }

      delete[] work;
      delete[] a_data;
      delete[] s_data;
      if (u_data) delete[] u_data;
      if (v_data) delete[] v_data;
    }


    /** Uses LAPACK's DSPSVX to solve symmetric non-definite linear system AX = B, where B is a RefSCVector
    */
    void lapack_linsolv_symmnondef(const RefSymmSCMatrix& A,
                                   RefSCVector& X,
                                   const RefSCVector& B)
    {
      int n = A.n();
      if (n != B.n())
        throw std::runtime_error("lapack_linsolv_symmnondef() -- dimensions of A and B do not match");
      if (n != X.n())
        throw std::runtime_error("lapack_linsolv_symmnondef() -- dimensions of A and X do not match");

      // convert A to packed upper triangular form
      int ntri = n*(n+1)/2;
      double* AP = new double[ntri];
      A.convert(AP);
      double* AFP = new double[ntri];
      int* ipiv = new int[n];
      double* BB = new double[n];
      B.convert(BB);
      double* XX = new double[n];

      char fact = 'N';
      char uplo = 'U';
      int nrhs = 1;
      double rcond = 0.0;
      double* ferr = new double[1];
      double* berr = new double[1];
      double* work = new double[3*n];
      int* iwork = new int[n];
      int info = 0;

      F77_DSPSVX(&fact, &uplo, &n, &nrhs, AP, AFP, ipiv, BB, &n, XX, &n, &rcond, ferr, berr, work, iwork, &info);

      if (info) {

        delete[] ferr;
        delete[] berr;
        delete[] work;
        delete[] iwork;

        delete[] AP;
        delete[] AFP;
        delete[] BB;
        delete[] XX;

        if (info > 0 && info <= n)
          throw std::runtime_error("lapack_linsolv_symmnondef() -- matrix A has factors which are exactly zero");
        if (info == n + 1)
          throw std::runtime_error("lapack_linsolv_symmnondef() -- matrix A is singular");
        if (info < 0)
          throw std::runtime_error("lapack_linsolv_symmnondef() -- one of the arguments to F77_DSPSVX is invalid");
      }
      else {
        X.assign(XX);

        delete[] ferr;
        delete[] berr;
        delete[] work;
        delete[] iwork;

        delete[] AP;
        delete[] AFP;
        delete[] BB;
        delete[] XX;
      }
    }

    /** Uses LAPACK's DSPSVX to solve symmetric non-definite linear system AX = B, where B is a RefSCMatrix.
        Dimensions of X and B must match.
    */
    void lapack_linsolv_symmnondef(const RefSymmSCMatrix& A,
                                   RefSCMatrix& X,
                                   const RefSCMatrix& B)
    {
      int n = A.n();
      int nrhs = B.ncol();
      if (n != B.nrow())
        throw std::runtime_error("lapack_linsolv_symmnondef() -- dimensions of A and B do not match");
      if (n != X.nrow())
        throw std::runtime_error("lapack_linsolv_symmnondef() -- dimensions of A and X do not match");
      if (X.ncol() != nrhs)
        throw std::runtime_error("lapack_linsolv_symmnondef() -- dimensions of B and X do not match");

      // convert A to packed upper triangular form
      int ntri = n*(n+1)/2;
      double* AP = new double[ntri];
      A.convert(AP);
      double* AFP = new double[ntri];
      int* ipiv = new int[n];
      double* BB = new double[nrhs*n];
      B.t().convert(BB);
      double* XX = new double[nrhs*n];

      char fact = 'N';
      char uplo = 'U';
      double rcond = 0.0;
      double* ferr = new double[nrhs];
      double* berr = new double[nrhs];
      double* work = new double[3*n];
      int* iwork = new int[n];
      int info = 0;

      F77_DSPSVX(&fact, &uplo, &n, &nrhs, AP, AFP, ipiv, BB, &n, XX, &n, &rcond, ferr, berr, work, iwork, &info);

      if (info) {

        delete[] ferr;
        delete[] berr;
        delete[] work;
        delete[] iwork;

        delete[] AP;
        delete[] AFP;
        delete[] BB;
        delete[] XX;

        if (info > 0 && info <= n)
          throw std::runtime_error("lapack_linsolv_symmnondef() -- matrix A has factors which are exactly zero");
        if (info == n + 1)
          throw std::runtime_error("lapack_linsolv_symmnondef() -- matrix A is singular");
        if (info < 0)
          throw std::runtime_error("lapack_linsolv_symmnondef() -- one of the arguments to F77_DSPSVX is invalid");
      }
      else {
        RefSCMatrix Xt = X.t();
        Xt.assign(XX);
        X = Xt.t();
        Xt = 0;

        delete[] ferr;
        delete[] berr;
        delete[] work;
        delete[] iwork;

        delete[] AP;
        delete[] AFP;
        delete[] BB;
        delete[] XX;
      }
    }


    /** Same as above, except uses C-style arrays already.
        A is in packed upper-triangular form, nA is the dimension of A,
        X and B are matrices with ncolB rows and nA columns.
    */
    void lapack_linsolv_symmnondef(const double* AP,
                                   int nA,
                                   double* Xt,
                                   const double* Bt,
                                   int ncolB)
    {
      const int n = nA;
      const int nrhs = ncolB;

      const int ntri = n*(n+1)/2;
      // AFP is a scratch array
      double* AFP = new double[ntri];
      int* ipiv = new int[n];

      char fact = 'N';
      char uplo = 'U';
      double rcond = 0.0;
      double* ferr = new double[nrhs];
      double* berr = new double[nrhs];
      double* work = new double[3*n];
      int* iwork = new int[n];
      int info = 0;

      F77_DSPSVX(&fact, &uplo, &n, &nrhs, AP, AFP, ipiv, Bt, &n, Xt, &n, &rcond, ferr, berr, work, iwork, &info);

      delete[] ferr;
      delete[] berr;
      delete[] work;
      delete[] iwork;

      delete[] AFP;

      if (info) {

        if (info > 0 && info <= n)
          throw std::runtime_error("lapack_linsolv_symmnondef() -- matrix A has factors which are exactly zero");
        if (info == n + 1)
          throw std::runtime_error("lapack_linsolv_symmnondef() -- matrix A is singular");
        if (info < 0)
          throw std::runtime_error("lapack_linsolv_symmnondef() -- one of the arguments to F77_DSPSVX is invalid");
      }
    }

    /** Computed eigenvalues of A and returns how many are below threshold. Uses LAPACK's DSYEVD.
    */
    int lapack_evals_less_eps(const RefSymmSCMatrix& A, double threshold)
    {
      int n = A.n();
      // convert A to square form
      double* Asq = new double[n*n];
      double** Arows = new double*[n];
      int ioff = 0;
      for(int i=0; i<n; i++, ioff+=n)
        Arows[i] = Asq + ioff;
      A.convert(Arows);
      delete[] Arows;
      double* evals = new double[n];

      char need_evecs = 'N';
      char uplo = 'U';
      const int lwork = 2*n+1;
      double* work = new double[lwork];
      const int liwork = 1;
      int* iwork = new int[liwork];
      int info = 0;

      F77_DSYEVD(&need_evecs, &uplo, &n, Asq, &n, evals, work, &lwork, iwork, &liwork, &info);

      if (info) {

        delete[] work;
        delete[] iwork;
        delete[] Asq;
        delete[] evals;

        throw sc::ProgrammingError("lapack_evals_less_eps() -- call to LAPACK's DSYEVD failed",__FILE__,__LINE__);
      }

      int nevals_lessthan_threshold = 0;
      for(int i=0; i<n; i++) {
        if (evals[i] < threshold)
          nevals_lessthan_threshold++;
      }

      delete[] work;
      delete[] iwork;
      delete[] Asq;
      delete[] evals;

      return nevals_lessthan_threshold;
    }

    void
    lapack_invert_symmnondef(RefSymmSCMatrix& A, double condition_number_threshold)
    {
      const int n = A.dim().n();
      const int ntri = n*(n+1)/2;
      char uplo = 'U';
      std::vector<double> AF(ntri);
      std::vector<int> ipiv(n);
      int info;
      std::vector<double> work(2*n);

      // compute factorization of A
      lapack_dpf_symmnondef(A, &(AF[0]), &(ipiv[0]), condition_number_threshold);

      // compute the inverse
      F77_DSPTRI(&uplo, &n, &(AF[0]), &(ipiv[0]), &(work[0]), &info);
      if (info) {
        if (info < 0)
          throw std::runtime_error("lapack_invert_symmnondef() -- one of the arguments to F77_DSPTRI is invalid");
        if (info > 0)
          throw std::runtime_error("lapack_invert_symmnondef() -- matrix A has factors which are exactly zero");
        assert(false);  // unreachable
      }

      A.assign( &(AF[0]) );
    }

    void lapack_dpf_symmnondef(const RefSymmSCMatrix& A,
                               double* AF,
                               int* ipiv,
                               double condition_number_threshold) {

      const int n = A.dim().n();
      char uplo = 'U';
      int info;
      std::vector<double> work(2*n);
      std::vector<int> iwork(n);

      // compute the infinity-norm of A
      A.convert(AF);
      char norm = 'I';
      const double anorm = F77_DLANSP(&norm, &uplo, &n, AF, &(work[0]));

      // factorize A = U.D.Ut
      F77_DSPTRF(&uplo, &n, AF, ipiv, &info);
      if (info) {
        if (info < 0)
          throw std::runtime_error("lapack_dpf_symmnondef() -- one of the arguments to F77_DSPTRF is invalid");
        if (info > 0)
          throw std::runtime_error("lapack_dpf_symmnondef() -- matrix A has factors which are exactly zero");
        assert(false);  // unreachable
      }

      // estimate the condition number
      double rcond;
      F77_DSPCON(&uplo, &n, AF, ipiv, &anorm, &rcond, &(work[0]), &(iwork[0]), &info);
      if (info) {
        if (info < 0)
          throw std::runtime_error("lapack_dpf_symmnondef() -- one of the arguments to F77_DSPCON is invalid");
        assert(false);  // unreachable
      }
      // if the condition number is above the threshold or its inverse below the working precision, throw
      if (condition_number_threshold != 0.0) {
        if (condition_number_threshold < 0.0)
          ExEnv::out0() << indent << "condition number estimate in lapack_dpf_symmnondef() = " << 1.0/rcond << std::endl;
        else if (1.0/rcond > condition_number_threshold)
          ExEnv::err0() << indent << "WARNING: large condition number in lapack_dpf_symmnondef(): threshold = "
                        << condition_number_threshold << " actual = " << 1.0/rcond << std::endl;
        const char epsilon = 'E';
        if (rcond < F77_DLAMCH(&epsilon))
          ExEnv::err0() << indent << "WARNING: condition number in lapack_dpf_symmnondef() exceeds the working precision" << std::endl;
      }

    }

    void
    lapack_linsolv_dpf_symmnondef(const double* A,
                                  int nA,
                                  const double* AF,
                                  const int* ipiv,
                                  double* Xt,
                                  const double* Bt,
                                  int ncolB)
    {
      const int n = nA;
      const char uplo = 'U';
      int info;

      // copy B into X
      const char full = 'F';
      F77_DLACPY(&full, &n, &ncolB, Bt, &n, Xt, &n, &info);

      // solve the linear system
      F77_DSPTRS(&uplo, &n, &ncolB, AF, ipiv, Xt, &n, &info);

      // Use iterative refinement to improve the computed solutions and
      // compute error bounds and backward error estimates for them.
      std::vector<double> ferr(ncolB);
      std::vector<double> berr(ncolB);
      std::vector<double> work(3*n);
      std::vector<int> iwork(n);
      F77_DSPRFS(&uplo, &n, &ncolB, A, AF, ipiv, Bt, &n, Xt, &n, &(ferr[0]),
                 &(berr[0]), &(work[0]), &(iwork[0]), &info);

    }

  };

};

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
