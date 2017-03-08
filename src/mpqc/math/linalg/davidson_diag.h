//
// Created by ChongPeng on 2/20/17.
//

#ifndef SRC_MPQC_MATH_LINALG_DAVIDSON_DIAG_H_
#define SRC_MPQC_MATH_LINALG_DAVIDSON_DIAG_H_

#include "mpqc/math/external/eigen/eigen.h"
#include "mpqc/math/linalg/gram_schmidt.h"
#include "mpqc/util/misc/exenv.h"

#include <TiledArray/algebra/utils.h>
#include <mpqc/util/misc/assert.h>
#include <tiledarray.h>

namespace mpqc {

template <typename Tile, typename Policy>
inline void plus(TA::DistArray<Tile, Policy>& y,
                 const TA::DistArray<Tile, Policy>& x) {
  const std::string vars =
      TA::detail::dummy_annotation(y.trange().tiles_range().rank());
  y(vars) += x(vars);
}

/**
 * \brief Davidson Algorithm
 *
 * \tparam D array type
 *
 *
 */

template <typename D>
class DavidsonDiag {
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
    EigenPair(const element_type&& value, const result_type&& vector)
        : eigen_value(std::move(value)), eigen_vector(std::move(vector)) {}

    ~EigenPair() = default;

    // sort by eigen value
    bool operator<(const EigenPair& other) const {
      return eigen_value < other.eigen_value;
    }
  };

 public:
  /**
   *
   * @param n_roots number of lowest roots to solve
   * @param n_guess number of eigen vector per root at subspace collapse,
   * default is 2
   * @param max_n_guess max number of guess vector per root, default is 3
   * @param symmetric if matrix is symmetric
   */
  DavidsonDiag(unsigned int n_roots, bool symmetric = true,
               unsigned int n_guess = 2, unsigned int max_n_guess = 3)
      : n_roots_(n_roots),
        symmetric_(symmetric),
        n_guess_(n_guess),
        max_n_guess_(max_n_guess),
        eigen_vector_() {}

  ~DavidsonDiag() { eigen_vector_.clear(); }
  /**
   *
   * @tparam Pred preconditioner object, which has void Pred(const element_type
   * & e,
   * D& residual) to update residual
   *
   * @param HB product with A and guess vector
   * @param B  guess vector
   * @param pred preconditioner
   * @return B updated guess vector
   */
  template <typename Pred>
  EigenVector<element_type> extrapolate(value_type& HB, value_type& B,
                                        const Pred& pred) {
    TA_ASSERT(HB.size() == B.size());
    // size of vector
    const auto n_v = B.size();

    // subspace
    // dot_product will return a replicated Eigen Matrix
    RowMatrix<element_type> G(n_v, n_v);
    for (auto i = 0; i < n_v; ++i) {
      for (auto j = 0; j < n_v; ++j) {
        G(i, j) = dot_product(B[j], HB[i]);
      }
    }

    // do eigen solve locally
    result_type E(n_roots_);
    RowMatrix<element_type> C(n_v, n_roots_);

    // symmetric matrix
    if (symmetric_) {
      // this return eigenvalue and eigenvector
      Eigen::SelfAdjointEigenSolver<RowMatrix<element_type>> es(G);

      RowMatrix<element_type> v = es.eigenvectors();
      EigenVector<element_type> e = es.eigenvalues();

      if (es.info() != Eigen::Success) {
        throw AlgorithmException("Eigen::SelfAdjointEigenSolver Failed!\n",
                                 __FILE__, __LINE__);
      }

      //        std::cout << es.eigenvalues() << std::endl;

      E = e.segment(0, n_roots_);
      C = v.leftCols(n_roots_);
    }
    // non-symmetric matrix
    else {
      // unitary transform to upper triangular matrix
      Eigen::RealSchur<RowMatrix<element_type>> rs(G);

      RowMatrix<element_type> T = rs.matrixT();
      RowMatrix<element_type> U = rs.matrixU();

      if (rs.info() != Eigen::Success) {
        throw AlgorithmException("Eigen::RealSchur Failed!\n", __FILE__,
                                 __LINE__);
      }

      // do eigen solve on T
      Eigen::EigenSolver<RowMatrix<element_type>> es(T);

      // sort eigen values
      std::vector<EigenPair> eg;
      {
        RowMatrix<element_type> v = es.eigenvectors().real();
        EigenVector<element_type> e = es.eigenvalues().real();

        if (rs.info() != Eigen::Success) {
          throw AlgorithmException("Eigen::EigenSolver Failed!\n", __FILE__,
                                   __LINE__);
        }

        for (auto i = 0; i < n_v; ++i) {
          eg.emplace_back(e[i], v.col(i));
        }

        std::sort(eg.begin(), eg.end());
      }

      // obtain final eigen value and eigen vector
      for (auto i = 0; i < n_roots_; ++i) {
        E[i] = eg[i].eigen_value;
        C.col(i) = U * eg[i].eigen_vector;
      }
    }

    // compute eigen_vector at current iteration and store it
    // X(i) = B(i)*C(i)
    value_type X(n_roots_);
    for (auto i = 0; i < n_roots_; ++i) {
      X[i] = copy(B[i]);
      zero(X[i]);
      for (auto j = 0; j < n_v; ++j) {
        axpy(X[i], C(j, i), B[j]);
      }
    }

    // check the size, if exceed n_guess, pop oldest
    if (eigen_vector_.size() == n_guess_) {
      eigen_vector_.pop_front();
    }
    eigen_vector_.push_back(X);

    // compute residual
    value_type residual(n_roots_);
    for (auto i = 0; i < n_roots_; ++i) {
      // initialize residual as 0
      residual[i] = copy(B[i]);
      zero(residual[i]);
      const auto e_i = -E[i];

      for (auto j = 0; j < n_v; ++j) {
        D tmp = copy(residual[i]);
        zero(tmp);
        axpy(tmp, e_i, B[j]);
        plus(tmp, HB[j]);
        scale(tmp, C(j, i));
        plus(residual[i], tmp);
      }
    }

    // precondition
    for (auto i = 0; i < n_roots_; i++) {
      pred(E[i], residual[i]);
    }

    B.insert(B.end(), residual.begin(), residual.end());

    // subspace collapse
    // Journal of Computational Chemistry, 11(10), 1164–1168.
    // https://doi.org/10.1002/jcc.540111008
    if(B.size() > n_roots_*max_n_guess_){
      B.clear();
      B.insert(B.end(), residual.begin(), residual.end());
      for(auto& vector: eigen_vector_){
        B.insert(B.end(), vector.begin(), vector.end());
      }
      // orthognolize all vectors
      gram_schmidt(B);
      // call it second times
      gram_schmidt(B);
    }
    else{
      // orthognolize new residual with original B
      gram_schmidt(B, n_v);
      // call it twice
      gram_schmidt(B, n_v);
    }


#ifndef NDEBUG
    const auto k = B.size();
    const auto tolerance =
        std::numeric_limits<typename D::element_type>::epsilon() * 100;
    for (auto i = 0; i < k; ++i) {
      for (auto j = i; j < k; ++j) {
        const auto test = dot_product(B[i], B[j]);
        //        std::cout << "i= " << i << " j= " << j << " dot= " << test <<
        //        std::endl;
        if (i == j) {
          TA_ASSERT(test - 1.0 < tolerance);
        } else {
          TA_ASSERT(test < tolerance);
        }
      }
    }
#endif

    return E.segment(0, n_roots_);
  }

 private:
  unsigned int n_roots_;
  bool symmetric_;
  unsigned int n_guess_;
  unsigned int max_n_guess_;
  std::deque<value_type> eigen_vector_;
};

}  // namespace mpqc

#endif  // SRC_MPQC_MATH_LINALG_DAVIDSON_DIAG_H_
