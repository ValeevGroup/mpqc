
//
// orthog.h -- orthogonalize the basis set
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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

#ifndef _chemistry_qc_basis_orthog_h
#define _chemistry_qc_basis_orthog_h

#include <util/state/state.h>
#include <math/scmat/matrix.h>

namespace sc {

/// This class computes the orthogonalizing transform for a basis set.
class OverlapOrthog: virtual public SavableState {
  public:

    /// An enum for the types of orthogonalization.
    enum OrthogMethod { Symmetric=1, Canonical=2, GramSchmidt=3 };
    /// default is to use Symmetric orthogonalization
    static OrthogMethod default_orthog_method() { return Symmetric; }
    /// default orthog threshold is 1e-8
    static double default_lindep_tol() { return 1e-8; }

  private:
    int debug_;

    RefSCDimension dim_;
    RefSCDimension orthog_dim_;

    // The tolerance for linearly independent basis functions.
    // The intepretation depends on the orthogonalization method.
    double lindep_tol_;
    // The number of linearly dependent functions
    int nlindep_;
    // The orthogonalization method
    OrthogMethod orthog_method_;
    // The orthogonalization matrices
    RefSCMatrix orthog_trans_;
    RefSCMatrix orthog_trans_inverse_;
    // The maximum and minimum residuals from the orthogonalization
    // procedure.  The interpretation depends on the method used.
    // For symmetry and canonical, these are the min and max overlap
    // eigenvalues.  These are the residuals for the basis functions
    // that actually end up being used.
    double min_orthog_res_;
    double max_orthog_res_;

    void compute_overlap_eig(RefSCMatrix& overlap_eigvec,
                             RefDiagSCMatrix& overlap_isqrt_eigval,
                             RefDiagSCMatrix& overlap_sqrt_eigval);
    void compute_symmetric_orthog();
    void compute_canonical_orthog();
    void compute_gs_orthog();
    void compute_orthog_trans();

    // WARNING: after a SavableState save/restore, these two members will
    // be null.  There is really no need to store these anyway--should be
    // removed.
    RefSymmSCMatrix overlap_;
    Ref<SCMatrixKit> result_kit_; // this kit is used for the result matrices

  public:
    /** @param lindep_tolerance If non-negative, this specifies the degree of linear dependence
               tolerated in the result. The precise definition is method-dependent. For symmetric/canonical
               methods this specifies the condition number of the truncated metric. If negative and
               method is symmetric/canonical,
               std::floor(-lindep_tolerance) specifies the number of linearly-dependent basis functions
               to be removed.
    */
    OverlapOrthog(OrthogMethod method,
                  const RefSymmSCMatrix &overlap,
                  const Ref<SCMatrixKit> &result_kit,
                  double lindep_tolerance,
                  int debug = 0);

    OverlapOrthog(StateIn&);

    virtual ~OverlapOrthog();

    void save_data_state(StateOut&);

    void reinit(OrthogMethod method,
                const RefSymmSCMatrix &overlap,
                const Ref<SCMatrixKit> &result_kit,
                double lindep_tolerance,
                int debug = 0);

    double min_orthog_res() const { return min_orthog_res_; }
    double max_orthog_res() const { return max_orthog_res_; }

    Ref<OverlapOrthog> copy() const;

    /// Returns the orthogonalization method
    OrthogMethod orthog_method() const { return orthog_method_; }

    /// Returns the tolerance for linear dependencies.
    double lindep_tol() const { return lindep_tol_; }

    /** Returns a matrix which does the requested transform from a basis to
        an orthogonal basis. The row dimension is the orthogonal basis
        dimension and the column dimension is the given basis dimension.

        This matrix can be used to transform operators from the given basis to the orthogonal basis
        and to transform functions from the orthogonal basis to the given basis.
        An operator in the orthogonal basis is obtained by \f$ X O_bas
        X^T\f$, where \f$X\f$ is the return value of this function and \f$ O_bas \f$
        is the operator in the given basis.
        A function in the given basis is obtained by \f$ X^T F_obas \f$, where
        \f$ F_obas \f$ is the function in the orthogonal basis.
        */
    RefSCMatrix basis_to_orthog_basis();

    /** Returns the inverse of the transformation returned by
        basis_to_orthog_basis() . The column dimension is the orthogonal basis
        dimension and the row dimension is the given basis dimension.

        This matrix can be used to transform operators from the orthogonal basis
        to the given basis and to transform functions from the orthogonal basis to the given basis.
        An operator in the given basis is obtained by \f$X O_obas
        X^T\f$, where \f$X\f$ is the return value of this function and \f$ O_obas \f$
        is the operator in the orthogonal basis.
        A function in the orthogonal basis is obtained by \f$ X^T F_bas \f$, where
        \f$ F_bas \f$ is the function in the given basis.
     */
    RefSCMatrix basis_to_orthog_basis_inverse();

    /** Return an $S^{-1}$. */
    RefSymmSCMatrix overlap_inverse();

    RefSCDimension dim();
    RefSCDimension orthog_dim();

    /** Returns the number of linearly dependent functions eliminated from
        the orthogonal basis.
      */
    int nlindep();
};

}

#endif
