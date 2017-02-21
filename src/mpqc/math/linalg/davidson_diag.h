//
// Created by ChongPeng on 2/20/17.
//

#ifndef MPQC_DAVIDSON_DIAG_H
#define MPQC_DAVIDSON_DIAG_H

#include "mpqc/math/external/eigen/eigen.h"
#include "mpqc/util/misc/exenv.h"

#include <TiledArray/algebra/utils.h>

namespace mpqc{


/**
 * \brief Davidson Algorithm
 *
 * \tparam D array type
 * \tparam F type that evaluates the LHS, which call F::operator()()
 *
 *
 */

template <typename D, typename F>
struct SymmDavidsonDiag {
  using value_type = typename D::element_type;
  using result_type = RowVector<value_type>;

  RowVector<value_type> operator() (F& H, D& B, std::size_t n_roots, value_type convergence, std::size_t max_iter){

    result_type result = result_type::Zero(n_roots);

    for(std::size_t iter = 0; iter < max_iter; ++iter){

      // compute H*B
      D HB = H(B);

      // subspace
      D G =  matrix_multiply(transpose(B), HB);

      result_type E;
      D V;
      std::tie(E, V) = eigen_solve(G);

      if(norm(E - result) < convergence){
        return E;
      }

      result = E;

      D BV = matrix_multiply(B, V);

      HB = matrix_multiply(HB, BV);

      EV = davidson_diag_operation1(BV, E);

      HB = HB - EV;

      HB = davidson_diag_operation2(HB, E, H.diagonal());

      // update B
      B = orthognolize(B,HB);

    }

    // print warning
    ExEnv::out0() << "\nWarning! Davidson Diagonalization Exceed Max Iteration! \n ";
    return result;

  }


};

} // namespace mpqc


#endif //MPQC_DAVIDSON_DIAG_H
