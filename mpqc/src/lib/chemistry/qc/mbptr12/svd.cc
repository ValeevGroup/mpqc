//
// svd.cc
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

#ifdef __GNUG__
#pragma implementation
#endif

#include <cmath>
#include <stdexcept>
#include <util/misc/formio.h>
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
        V      U^t
        Sigma  Sigma^t
        U      Vt^t

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
        V = V.t();
      }
      
      delete[] work;
      delete[] a_data;
      delete[] s_data;
      if (u_data) delete[] u_data;
      if (v_data) delete[] v_data;
    }

  };
};

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
