//
// Created by ChongPeng on 2/20/17.
//

#ifndef MPQC_DAVIDSON_DIAG_H
#define MPQC_DAVIDSON_DIAG_H

#include "mpqc/chemistry/qc/lcao/expression/trange1_engine.h"
#include "mpqc/math/external/eigen/eigen.h"
#include "mpqc/math/tensor/clr/array_to_eigen.h"
#include "mpqc/util/misc/exenv.h"

#include <TiledArray/algebra/utils.h>
#include <mpqc/util/misc/assert.h>
#include <tiledarray.h>

namespace mpqc {

template <typename Tile, typename Policy>
inline void plus(TA::DistArray<Tile,Policy>& y,
                 const TA::DistArray<Tile,Policy>& x) {
  const std::string vars = TA::detail::dummy_annotation(y.trange().tiles_range().rank());
  y(vars) = y(vars) + x(vars);
}

/**
 * \brief Davidson Algorithm
 *
 * \tparam D array type
 *
 *
 */

template <typename D>
class SymmDavidsonDiag {
 public:
  using element_type = typename D::element_type;
  using result_type = EigenVector<element_type>;
  using value_type = std::vector<D>;

 private:
  struct EigenPair {
    element_type eigen_value;
    result_type eigen_vector;

    /// constructor
    EigenPair(const element_type& value, const result_type& vector)
        : eigen_value(value), eigen_vector(vector) {}

    /// move constructor
//    EigenPair(const element_type&& value, const result_type&& vector)
//        : eigen_value(std::move(value)), eigen_vector(std::move(vector)) {}

    ~EigenPair() = default;

    // sort by eigen value
    bool operator<(const EigenPair& other) const {
      return eigen_value < other.eigen_value;
    }
  };

 public:
  SymmDavidsonDiag(unsigned int n_roots, unsigned int n_guess)
      : n_roots_(n_roots), n_guess_(n_guess) {}

  EigenVector<element_type> extrapolate(value_type& HB, value_type& B) {
    TA_ASSERT(HB.size() == B.size());
    // size of vector
    const auto n_v = B.size();

    // subspace
    // dot_product will return a replicated Eigen Matrix
    RowMatrix<element_type> G(n_v, n_v);
    for(auto i = 0; i < n_v; ++i){
      for(auto j = 0; j < n_v; ++j){
        G(i,j) = dot_product(B[j], HB[i]);
      }
    }

    // do eigen solve locally
    result_type E(n_guess_);
    RowMatrix<element_type> C(n_v, n_guess_);
    {
      // this return eigenvalue and eigenvector with size n_guess
      Eigen::EigenSolver<RowMatrix<element_type>> es(G);

      // sort eigenvalue and eigenvector
      std::vector<EigenPair> eg;
      {
        RowMatrix<element_type> v = es.eigenvectors().real();
        EigenVector<element_type> e = es.eigenvalues().real();

        for (auto i = 0; i < n_guess_; ++i) {
          eg.emplace_back(e[i], v.col(i));
        }

        std::sort(eg.begin(), eg.end());
      }

      // obtain eigenvalue and eigenvector

      for (auto i = 0; i < n_guess_; ++i) {
        E[i] = eg[i].eigen_value;
        C.col(i) = eg[i].eigen_vector;
      }

    }

    // compute residual
    value_type residual(n_guess_);
    for(auto i = 0; i < n_guess_; ++i){

      // initialize redidual as 0
      residual[i] = copy(B[i]);
      zero(residual[i]);
      const auto e_i = -E[i];
      for(auto j = 0; j < n_v; ++j){

        D tmp = copy(residual[i]);
        zero(tmp);
        axpy(tmp,e_i,B[i]);
        plus(tmp,HB[i]);
        scale(tmp, C(j,i));
        plus(residual[i], tmp);

      }
    }

    // orthognolize residual with B
    B.insert(B.end(), residual.begin(), residual.end());

    return E.segment(0, n_roots_);
  }

 private:
  unsigned int n_roots_;
  unsigned int n_guess_;
};

}  // namespace mpqc

#endif  // MPQC_DAVIDSON_DIAG_H
