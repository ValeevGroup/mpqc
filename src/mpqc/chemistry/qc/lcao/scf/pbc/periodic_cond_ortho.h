#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_COND_ORTHO_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_COND_ORTHO_H_

#include "mpqc/chemistry/qc/lcao/factory/periodic_ao_factory.h"
#include "mpqc/math/external/eigen/util.h"

namespace mpqc {
namespace utility {

/**
 * \brief gensqrtinv
 *
 * \param S is complex overlap matrix for a specific k vector
 * \param S_condition_num_thresh is maximum condition number allowed
 * \return
 */
template <typename TArray = TA::DistArray<TA::TensorZ, TA::SparsePolicy>>
MatrixZ gensqrtinv(const TArray &S, bool symmetric, double max_condition_num,
                   int64_t k, bool print = false) {
  auto& world = S.world();

  auto S_eig = math::array_to_eigen(S);
  assert(S_eig.rows() == S_eig.cols());

  MatrixZ X;
  if (world.rank() == 0) {

    // compute generalized orthogonalizer X such that Xt.S.X = I
    // if symmetric, X = U.s_sqrtinv.Ut
    // if canonical, X = U.s_sqrtinv
    // where s and U are eigenvalues and eigenvectors of S
    Eigen::SelfAdjointEigenSolver<MatrixZ> comp_eig_solver(S_eig);
    auto U = comp_eig_solver.eigenvectors();
    auto s = comp_eig_solver.eigenvalues();
    //  ::mpqc::detail::sort_eigen(s, U);

    auto s_max = s.maxCoeff();
    auto s_min = s.minCoeff();

    auto S_condition_num =
        std::min(s_max / std::max(s_min, std::numeric_limits<double>::min()),
                 1.0 / std::numeric_limits<double>::epsilon());

    auto threshold = s_max / max_condition_num;

    int64_t s_rows = s.rows();
    int64_t s_cond = 0;
    for (int64_t i = s_rows - 1; i >= 0; --i) {
      if (s(i) >= threshold) {
        ++s_cond;
      } else
        i = 0;
    }

    auto sigma = s.bottomRows(s_cond);
    auto result_condition_num = sigma.maxCoeff() / sigma.minCoeff();
    auto sigma_invsqrt_real =
        Eigen::MatrixXd(sigma.array().sqrt().inverse().matrix().asDiagonal());
    auto sigma_invsqrt = sigma_invsqrt_real.cast<std::complex<double>>();

    // make canonical X
    auto U_cond = U.block(0, s_rows - s_cond, s_rows, s_cond);
    X = U_cond * sigma_invsqrt;
    // make symmetric X
    if (symmetric) {
      X = X * U_cond.transpose().conjugate();
    }

    auto nbf_omitted = s_rows - s_cond;
    if (nbf_omitted < 0) throw "Error: dropping negative number of functions!";

    if (print) {
      ExEnv::out0() << "\n\toverlap condition number = " << S_condition_num
                    << " at k = " << k;
    }

    if (nbf_omitted > 0) {
      auto XtS = X.transpose().conjugate() * S_eig;
      auto should_be_I = XtS * X;

      auto I_real =
          Eigen::MatrixXd::Identity(should_be_I.rows(), should_be_I.cols());
      auto I_comp = I_real.template cast<std::complex<double>>();
      auto should_be_zero = (should_be_I - I_comp).norm();

      ExEnv::out0() << " (dropped " << nbf_omitted << " "
                    << (nbf_omitted > 1 ? "fns" : "fn") << " to reduce to "
                    << result_condition_num << ")" << std::endl;
      ExEnv::out0() << "\t\t||Xt*S*X - I||_2 = " << should_be_zero
                    << " (should be zero)" << std::endl;
    }
  }  // rank == 0
  world.gop.broadcast_serializable(X, 0);

  return X;
}

/**
 * \brief conditioned_orthogonalizer
 *
 * \param overlap is rectangular overlap matrix in k_space
 * \param k_size is total # of k points
 * \param max_condition_num is maximum condition number allowed
 * \return conditioned orthogonalizer
 */
template <typename TArray = TA::DistArray<TA::TensorZ, TA::SparsePolicy>>
std::vector<MatrixZ> conditioned_orthogonalizer(
    const TArray overlap, int64_t k_size, double max_condition_num = 1.0e8,
    int64_t print_max_item = 100) {
  std::vector<MatrixZ> X;
  X.resize(k_size);

  auto tr0 = overlap.trange().data()[0];
  auto tr1 = overlap.trange().data()[1];

  assert(tr1.extent() == (tr0.extent() * k_size));

  auto n_tr_nu = overlap.trange().data().front().tiles_range().second;

  bool print = (k_size < print_max_item) ? true : false;
  for (auto k = 0; k < k_size; ++k) {
    std::vector<std::size_t> S_low{0, k * n_tr_nu};
    std::vector<std::size_t> S_up{n_tr_nu, (k + 1) * n_tr_nu};

    TArray S;
    S("mu, nu") = overlap("mu, nu").block(S_low, S_up);
    X[k] = gensqrtinv(S, false, max_condition_num, k, print);
  }

  ExEnv::out0() << std::endl;
  return X;
}

}  // namespace utility
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_COND_ORTHO_H_
