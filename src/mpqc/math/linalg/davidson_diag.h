//
// Created by ChongPeng on 2/20/17.
//

#ifndef SRC_MPQC_MATH_LINALG_DAVIDSON_DIAG_H_
#define SRC_MPQC_MATH_LINALG_DAVIDSON_DIAG_H_

#include <TiledArray/algebra/utils.h>
#include <tiledarray.h>

#include "mpqc/math/external/eigen/eigen.h"
#include "mpqc/math/linalg/gram_schmidt.h"
#include "mpqc/util/misc/assert.h"
#include "mpqc/util/misc/exception.h"
#include "mpqc/util/misc/exenv.h"
#include "mpqc/util/misc/print.h"
#include "mpqc/util/misc/time.h"

namespace mpqc {

template <typename D>
class DavidsonDiagPred {
 public:
  virtual void operator()(const EigenVector<typename D::element_type>& e,
                          std::vector<D>& guess) const = 0;

  virtual typename D::element_type norm(const D& d) const {
    return norm2(d) / size(d);
  }

  virtual ~DavidsonDiagPred() = default;
};

// clang-format off
/**
 * \brief Davidson Algorithm
 *
 * Solves eigen value problem <tt> Hx = ex </tt> for the lowest n eigen value and eigen vector
 *
 * it starts with orthogonal guess vector B {b1, b2, ... bn}
 * the eigen vector is a linear combination of B
 * extrapolate() will update the vector B and store new x
 *
 * \tparam D
 * array type
 * \c D::element_type must be defined and \c D must provide the following stand-alone functions:
 * - `D copy(const D&)`
 * - `element_type dot_product(const D& a, const D& b)`
 * - `void scale(D& y , element_type a)`
 * - `void axpy(D&y , element_tye a, const D& z)`
 * - `void zero(D& x)`
 * - `std::size_t size(const D& x)`
 * - `element_type norm2(const D& x)`
 *
 */
// clang-format on
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

   * @param symmetric if matrix is symmetric

   * @param n_guess number of eigen vector per root at subspace collapse,
   * default is 2

   * @param max_n_guess max number of guess vector per root, default is 4
   *
   * @param vector_threshold threshold for the norm of new guess vector in gram
   schmidt orthonormalization
   */
  DavidsonDiag(unsigned int n_roots, bool symmetric = true,
               unsigned int n_guess = 2, unsigned int max_n_guess = 4,
               double vector_threshold = 1.0e-5)
      : n_roots_(n_roots),
        symmetric_(symmetric),
        n_guess_(n_guess),
        max_n_guess_(max_n_guess),
        vector_threshold_(vector_threshold),
        eigen_vector_(),
        HB_(),
        B_(),
        subspace_() {}

  virtual ~DavidsonDiag() {
    eigen_vector_.clear();
    HB_.clear();
    B_.clear();
    subspace_.resize(0, 0);
  }

  /**
   *
   * @tparam Operator  operator that computes the product of H*B
   *
   * @param guess initial guess vector
   * @param op    op(B) should compute HB
   * @param pred  preconditioner, which inherit from DavidsonDiagPred
   * @param convergence   convergence threshold
   * @param max_iter  max number of iteration allowd
   * @return
   */
  template <typename Operator>
  EigenVector<element_type> solve(value_type& guess, const Operator& op,
                                  const DavidsonDiagPred<D>* const pred,
                                  double convergence, std::size_t max_iter) {
    double norm_e = 1.0;
    double norm_r = 1.0;
    std::size_t iter = 0;
    auto& world = TA::get_default_world();

    EigenVector<element_type> eig = EigenVector<element_type>::Zero(n_roots_);

    while (iter < max_iter && (norm_r > convergence || norm_e > convergence)) {
      auto time0 = mpqc::fenced_now(world);

      std::size_t dim = guess.size();
      //    ExEnv::out0() << "vector dimension: " << dim << std::endl;

      // compute product of H with guess vector
      value_type HC = op(guess);

      auto time1 = mpqc::fenced_now(world);
      EigenVector<element_type> eig_new, norms;
      std::tie(eig_new, norms) = extrapolate(HC, guess, pred);
      auto time2 = mpqc::fenced_now(world);

      EigenVector<element_type> delta_e = (eig - eig_new);
      delta_e = delta_e.cwiseAbs();
      norm_e =
          *std::max_element(delta_e.data(), delta_e.data() + delta_e.size());
      norm_r = *std::max_element(norms.data(), norms.data() + norms.size());

      util::print_excitation_energy_iteration(
          iter, delta_e, norms, eig_new, mpqc::duration_in_s(time0, time1),
          mpqc::duration_in_s(time1, time2));

      eig = eig_new;
      iter++;

    }  // end of while loop

    if (iter == max_iter) {
      throw MaxIterExceeded("Davidson Diagonalization Exceeded Max Iteration",
                            __FILE__, __LINE__, max_iter, "DavidsonDiag");
    }

    return eig;
  };

  /// @return return current eigen vector in Davidson
  virtual value_type& eigen_vector() { return eigen_vector_.back(); }

  // clang-format off
  /**
   *
   * @param HB product with A and guess vector
   * @param B  guess vector
   * @param pred preconditioner, which inherit from DavidsonDiagPred
   *
   * @return B updated guess vector
   * @return updated eigen values, norm of residual
   */
  // clang-format on
  std::tuple<EigenVector<element_type>, EigenVector<element_type>> extrapolate(
      value_type& HB, value_type& B, const DavidsonDiagPred<D>* const pred) {
    TA_ASSERT(HB.size() == B.size());
    // size of new vector
    const auto n_b = B.size();
    // size of original subspace
    const auto n_s = subspace_.cols();

    B_.insert(B_.end(), B.begin(), B.end());
    B.clear();

    HB_.insert(HB_.end(), HB.begin(), HB.end());
    HB.clear();
    // size of new subspace
    const auto n_v = B_.size();

    // compute new subspace
    // G will be replicated Eigen Matrix
    {
      RowMatrix<element_type> G = RowMatrix<element_type>::Zero(n_v, n_v);
      // reuse stored subspace
      G.block(0, 0, n_s, n_s) << subspace_;
      // initialize new value
      if (symmetric_) {
        for (std::size_t i = 0; i < n_b; ++i) {
          const auto ii = i + n_s;
          for (std::size_t j = 0; j <= ii; ++j) {
            G(ii, j) = dot_product(B_[ii], HB_[j]);
            if (ii != j) {
              G(j, ii) = G(ii, j);
            }
          }
        }
      } else {
        for (std::size_t i = 0; i < n_b; ++i) {
          const auto ii = i + n_s;
          for (std::size_t j = 0; j <= ii; ++j) {
            G(ii, j) = dot_product(B_[ii], HB_[j]);
            if (ii != j) {
              G(j, ii) = dot_product(B_[j], HB_[ii]);
            }
          }
        }
      }
      subspace_ = G;
    }

    //    std::cout << "G: " << std::endl;
    //    std::cout << G << std::endl;

    deflation(subspace_);

    // do eigen solve locally
    result_type E(n_roots_);
    RowMatrix<element_type> C(n_v, n_roots_);

    // symmetric matrix
    if (symmetric_) {
      // this return eigenvalue and eigenvector
      Eigen::SelfAdjointEigenSolver<RowMatrix<element_type>> es(subspace_);

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
      // do eigen solve on T
      Eigen::EigenSolver<RowMatrix<element_type>> es(subspace_);

      // sort eigen values
      std::vector<EigenPair> eg;
      {
        RowMatrix<element_type> v = es.eigenvectors().real();
        EigenVector<element_type> e = es.eigenvalues().real();

        if (es.info() != Eigen::Success) {
          throw AlgorithmException("Eigen::EigenSolver Failed!\n", __FILE__,
                                   __LINE__);
        }

        for (std::size_t i = 0; i < n_v; ++i) {
          eg.emplace_back(e[i], v.col(i));
        }

        std::sort(eg.begin(), eg.end());
      }

      // obtain final eigen value and eigen vector
      for (std::size_t i = 0; i < n_roots_; ++i) {
        E[i] = eg[i].eigen_value;
        C.col(i) = eg[i].eigen_vector;
      }

      // orthonormalize C
      //      RowMatrix<element_type> Q = C;
      //      Eigen::ColPivHouseholderQR<RowMatrix<element_type>> qr(Q);
      //      Eigen::ColPivHouseholderQR<RowMatrix<element_type>> qr(C);
      //      C = qr.householderQ();

      //      const auto tolerance =
      //          std::numeric_limits<typename D::element_type>::epsilon() *
      //          100;
      //      for (auto i = 0; i < n_roots_; ++i) {
      //        for (auto j = i; j < n_roots_; ++j) {
      //          const auto test = C.col(i).dot(C.col(j));
      //          if (i == j) {
      //            TA_ASSERT(test - 1.0 < tolerance);
      //          } else {
      //            TA_ASSERT(test < tolerance);
      //          }
      //          std::cout << "i= " << i << " j= " << j << " dot= " << test
      //                    << std::endl;
      //        }
      //      }
    }

    // compute eigen_vector at current iteration and store it
    // X(i) = B(i)*C(i)
    value_type X(n_roots_);
    for (std::size_t i = 0; i < n_roots_; ++i) {
      X[i] = copy(B_[i]);
      zero(X[i]);
      for (std::size_t j = 0; j < n_v; ++j) {
        axpy(X[i], C(j, i), B_[j]);
      }
    }

    // check the size, if exceed n_guess, pop oldest
    if (eigen_vector_.size() == n_guess_) {
      eigen_vector_.pop_front();
    }
    eigen_vector_.push_back(X);

    // compute residual
    // R(i) = (H - e(i)I)*B(i)*C(i)
    //      = (HB(i)*C(i) - e(i)*X(i)
    value_type residual(n_roots_);
    EigenVector<element_type> norms(n_roots_);
    for (std::size_t i = 0; i < n_roots_; ++i) {
      residual[i] = copy(X[i]);
      const auto e_i = -E[i];
      scale(residual[i], e_i);
      for (std::size_t j = 0; j < n_v; ++j) {
        axpy(residual[i], C(j, i), HB_[j]);
      }
      norms[i] = pred->norm(residual[i]);
    }

    // precondition
    // user should define preconditioner
    // usually it is D(i) = (e(i) - H_D)^-1 R(i)
    // where H_D is the diagonal element of H
    // but H_D can be approximated and computed on the fly
    pred->operator()(E, residual);

    // subspace collapse
    // restart with new vector and most recent eigen vector
    // Journal of Computational Chemistry, 11(10), 1164–1168.
    // https://doi.org/10.1002/jcc.540111008
    if (B_.size() > n_roots_ * (max_n_guess_ - 1)) {
      B_.clear();
      HB_.clear();
      subspace_.resize(0, 0);
      B.insert(B.end(), residual.begin(), residual.end());
      // TODO this generates too much eigen vectors
      // use all stored eigen vector from last n_guess interation
      for (auto& vector : eigen_vector_) {
        B.insert(B.end(), vector.begin(), vector.end());
      }
      // orthognolize all vectors
      gram_schmidt(B, vector_threshold_);
      // call it second times
      gram_schmidt(B, vector_threshold_);
    } else {
      // TODO better way to orthonormalize than double gram_schmidt
      // orthognolize new residual with original B
      gram_schmidt(B_, residual, vector_threshold_);
      // call it twice
      gram_schmidt(B_, residual, vector_threshold_);
      B = residual;

      //      for (std::size_t i = 0; i < n_roots_; i++) {
      //        const auto m = norm2(B[i]);
      //        std::cout << "norm: " << i << " " << m << std::endl;
      //        TA_ASSERT(m > 1.0e-3);
      //      }
    }

// test if orthonomalized
#ifndef NDEBUG
    const auto k = B.size();
    const auto tolerance =
        std::numeric_limits<typename D::element_type>::epsilon() * 100;
    for (std::size_t i = 0; i < k; ++i) {
      for (std::size_t j = i; j < k; ++j) {
        const auto test = dot_product(B[i], B[j]);
        //                std::cout << "i= " << i << " j= " << j << " dot= " <<
        //                test <<
        //                std::endl;
        if (i == j) {
          TA_ASSERT(test - 1.0 < tolerance);
        } else {
          TA_ASSERT(test < tolerance);
        }
      }
    }
#endif

    return std::make_tuple(E.segment(0, n_roots_), norms);
  }

 private:
  virtual void deflation(RowMatrix<element_type>& A) const {
    // do nothing
  }

 protected:
  unsigned int n_roots_;
  bool symmetric_;
  unsigned int n_guess_;
  unsigned int max_n_guess_;
  double vector_threshold_;
  std::deque<value_type> eigen_vector_;
  value_type HB_;
  value_type B_;
  RowMatrix<element_type> subspace_;
};

template <typename D>
class SingleStateDavidsonDiag : public DavidsonDiag<D> {
 public:
  using typename DavidsonDiag<D>::element_type;
  using typename DavidsonDiag<D>::value_type;
  using typename DavidsonDiag<D>::result_type;

  SingleStateDavidsonDiag(unsigned int n_roots, double shift,
                          bool symmetric = true, unsigned int n_guess = 2,
                          unsigned int max_n_guess = 4,
                          double vector_threshold = 1.0e-5)
      : DavidsonDiag<D>(n_roots, symmetric, n_guess, max_n_guess,
                     vector_threshold),
        total_roots_(n_roots),
        shift_(shift) {}

  /**
   * This is not a virtual function, it doesn't override DavidsonDiag::solve()
   *
   * @tparam Operator  operator that computes the product of H*B
   *
   * @param guess initial guess vector
   * @param op    op(B) should compute HB
   * @param pred  preconditioner, which inherit from DavidsonDiagPred
   * @param convergence   convergence threshold
   * @param max_iter  max number of iteration allowd
   * @return
   */
  template <typename Operator>
  EigenVector<element_type> solve(value_type& guess, const Operator& op,
                                  const DavidsonDiagPred<D>* const pred,
                                  double convergence, std::size_t max_iter) {
    // set roots in DavidsonDiag as 1, solve roots 1 at a time
    this->n_roots_ = 1;

    EigenVector<element_type> total_eig =
        EigenVector<element_type>::Zero(total_roots_);

    TA_ASSERT(guess.size() == total_roots_);

    for (std::size_t i = 0; i < total_roots_; i++) {
      ExEnv::out0() << "Start solved root: " << i << "\n";
      std::vector<D> guess_i = {guess[i]};

      double norm_e = 1.0;
      double norm_r = 1.0;
      std::size_t iter = 0;
      auto& world = TA::get_default_world();

      EigenVector<element_type> eig = EigenVector<element_type>::Zero(1);

      while (iter < max_iter &&
             (norm_r > convergence || norm_e > convergence)) {
        auto time0 = mpqc::fenced_now(world);

        // compute product of H with guess vector
        value_type HC = op(guess_i);

        auto time1 = mpqc::fenced_now(world);
        EigenVector<element_type> eig_new, norms;
        std::tie(eig_new, norms) = this->extrapolate(HC, guess_i, pred);
        auto time2 = mpqc::fenced_now(world);

        EigenVector<element_type> delta_e = (eig - eig_new);
        delta_e = delta_e.cwiseAbs();
        norm_e = delta_e[0];
        norm_r = norms[0];

        util::print_excitation_energy_iteration(
            iter, delta_e, norms, eig_new, mpqc::duration_in_s(time0, time1),
            mpqc::duration_in_s(time1, time2));

        eig = eig_new;
        iter++;

      }  // end of while loop

      if (iter == max_iter) {
        throw MaxIterExceeded("Davidson Diagonalization Exceeded Max Iteration",
                              __FILE__, __LINE__, max_iter,
                              "SingleStateDavidsonDiag");
      }

      // converged
      converged_eigen_vector_.push_back(this->eigen_vector()[0]);
      reset();

      total_eig[i] = eig[0];
      ExEnv::out0() << "\n";
      ExEnv::out0() << "Solved root: " << i << " value:" << eig[0] << "\n";
    }
  }

 private:
  void reset() {
    this->eigen_vector_.clear();
    this->HB_.clear();
    this->B_.clear();
    this->subspace_.resize(0, 0);
  }

  void deflation(RowMatrix<element_type>& A) const override {
//    std::size_t n = converged_eigen_vector_.size();

//    element_type shift = 0.0;

//    for (std::size_t i = 0; i < n; i++) {
//      shift +=
//          dot_product(converged_eigen_vector_[i], converged_eigen_vector_[i]);
//    }

//    std::cout << "shift: " << shift << "\n";

//    A = A.array() + shift_ * shift;
  }

  std::size_t total_roots_;
  element_type shift_;
  value_type converged_eigen_vector_;
};

}  // namespace mpqc

#endif  // SRC_MPQC_MATH_LINALG_DAVIDSON_DIAG_H_
