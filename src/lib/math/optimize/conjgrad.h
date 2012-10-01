//
// conjgrad.h
//
// Copyright (C) 2012 Edward Valeev
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
#pragma interface
#endif

#ifndef _mpqc_src_lib_math_optimize_conjgrad_h
#define _mpqc_src_lib_math_optimize_conjgrad_h

#include <assert.h>
#include <math/scmat/abstract.h>

namespace sc {

  /**
   * Solves linear system a(x) = b using conjugate gradient solver
   * where a is a linear function of x.
   *
   * @param a functor that evaluates LHS, will call as: a(x, result)
   * @param b RHS
   * @param x unknown
   * @param preconditioner
   * @param convergence_target
   * @return
   */
  template <typename F>
  double linsolv_conjugate_gradient(F& a,
                                    const RefSCMatrix& b,
                                    RefSCMatrix& x,
                                    const RefSCMatrix& preconditioner,
                                    double convergence_target = -1.0) {

    const size_t n = x.nrow() * x.ncol();
    assert(n == (preconditioner.nrow() * preconditioner.ncol()));

    // solution vector
    RefSCMatrix XX_i;
    // residual vector
    RefSCMatrix RR_i = b.clone();
    // preconditioned residual vector
    RefSCMatrix ZZ_i;
    // direction vector
    RefSCMatrix PP_i;
    RefSCMatrix APP_i = b.clone();

    // approximate the condition number as the ratio of the min and max elements of the preconditioner
    // assuming that preconditioner is the approximate inverse of A in Ax - b =0
    typedef SCElementFindExtremum<std::greater<double>, SCMatrixIterationRanges::AllElements> MinOp;
    typedef SCElementFindExtremum<std::less<double>, SCMatrixIterationRanges::AllElements> MaxOp;
    Ref<MinOp> findmin_op = new MinOp;
    Ref<MaxOp> findmax_op = new MaxOp;
    preconditioner.element_op(findmin_op);
    const double precond_min = findmin_op->result().at(0).value;
    const double precond_max = findmax_op->result().at(0).value;
    const double cond_number = precond_max / precond_min;
    if (convergence_target < 0.0) convergence_target = 1e-15 * cond_number; // estimate of how tightly the system can be converged

    bool converged = false;
    const unsigned int max_niter = 500;
    double rnorm2 = 0.0;

    // starting guess: x_0 = D^-1 . b
    Ref<SCElementOp2> multiply_op = new SCDestructiveElementProduct;
    XX_i = b.copy();
    XX_i.element_op(multiply_op, preconditioner);

    // r_0 = b - a(x)
    a(XX_i, RR_i);  // RR_i = a(XX_i)
    RR_i.scale(-1.0);
    RR_i.accumulate(b); // RR_i = b - a(XX_i)

    // z_0 = D^-1 . r_0
    ZZ_i = RR_i.copy();
    ZZ_i.element_op(multiply_op, preconditioner);

    // p_0 = z_0
    PP_i = ZZ_i.copy();

    Ref<SCElementScalarProduct> scalarprod_op = new SCElementScalarProduct;

    unsigned int iter = 0;
    while (not converged) {

      // alpha_i = (r_i . z_i) / (p_i . A . p_i)
      scalarprod_op->init();
      RR_i.element_op(scalarprod_op, ZZ_i);
      const double rz_norm2 = scalarprod_op->result();
      a(PP_i,APP_i);

      scalarprod_op->init();
      PP_i.element_op(scalarprod_op, APP_i);
      const double pAp_i = scalarprod_op->result();
      const double alpha_i = rz_norm2 / pAp_i;

      // x_i += alpha_i p_i
      Ref<SCElementDAXPY> daxpy_op = new SCElementDAXPY(alpha_i);
      XX_i.element_op(daxpy_op, PP_i);

      // r_i -= alpha_i Ap_i
      Ref<SCElementDAXPY> daxpy_op2 = new SCElementDAXPY(-alpha_i);
      RR_i.element_op(daxpy_op2, APP_i);

      Ref<SCElementKNorm> norm2_op = new SCElementKNorm(2);
      RR_i.element_op(norm2_op);
      const double r_ip1_norm = norm2_op->result();
      if (r_ip1_norm < convergence_target) {
        converged = true;
        rnorm2 = r_ip1_norm;
      }

      // z_i = D^-1 . r_i
      ZZ_i.assign(RR_i);
      ZZ_i.element_op(multiply_op, preconditioner);
      scalarprod_op->init();
      ZZ_i.element_op(scalarprod_op, RR_i);
      const double rz_ip1_norm2 = scalarprod_op->result();

      const double beta_i = rz_ip1_norm2 / rz_norm2;

      // p_i = z_i+1 + beta_i p_i
      // 1) scale p_i by beta_i
      // 2) add z_i+1 (i.e. current contents of z_i)
      PP_i.scale(beta_i);
      Ref<SCElementDAXPY> daxpy_op3 = new SCElementDAXPY( 1.0 );
      PP_i.element_op(daxpy_op3, ZZ_i);

      ++iter;
      //std::cout << "iter=" << iter << " dnorm=" << r_ip1_norm << std::endl;

      if (iter >= max_niter) {
        x.assign(XX_i);
        throw MaxIterExceeded("linsolv_conjugate_gradient", __FILE__, __LINE__, max_niter);
      }
    } // solver loop

    x.assign(XX_i);

    return rnorm2;
  }

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
