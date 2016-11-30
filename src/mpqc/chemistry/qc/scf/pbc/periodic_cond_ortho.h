#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_CONDITIONED_ORTHOGONALIZER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_CONDITIONED_ORTHOGONALIZER_H_

#include "mpqc/chemistry/qc/integrals/periodic_ao_factory.h"

namespace mpqc {
namespace scf {

/**
 * \brief gensqrtinv
 *
 * \param S is complex overlap matrix for a specific k vector
 * \param S_condition_num_thresh is maximum condition number allowed
 * \return
 */
template <typename TArray = TA::DistArray<TA::TensorZ, TA::SparsePolicy>>
Matrixc gensqrtinv(const TArray S, bool symmetric, double max_condition_num, int64_t k) {
  double S_condition_num;
  Matrixc X;
  auto &world = S.world();
  auto S_eig = array_ops::array_to_eigen(S);

  assert(S_eig.rows() == S_eig.cols());

  // compute generalized orthogonalizer X such that Xt.S.X = I
  // if symmetric, X = U.s_sqrtinv.Ut
  // if canonical, X = U.s_sqrtinv
  // where s and U are eigenvalues and eigenvectors of S
  Eigen::ComplexEigenSolver<Matrixc> comp_eig_solver(S_eig);
  auto U = comp_eig_solver.eigenvectors();
  auto s = comp_eig_solver.eigenvalues();
  integrals::detail::sort_eigen(s, U);

  auto s_real_max = s.real().maxCoeff();
  auto s_real_min = s.real().minCoeff();

  S_condition_num = std::min(
      s_real_max / std::max(s_real_min, std::numeric_limits<double>::min()),
      1.0 / std::numeric_limits<double>::epsilon());

  auto threshold = s_real_max / max_condition_num;

  int64_t s_rows = s.rows();
  int64_t s_cond = 0;
  for (int64_t i = s_rows - 1; i >= 0; --i) {
    if (s.real()(i) >= threshold) {
      ++s_cond;
    } else
      i = 0;
  }

  auto sigma = s.bottomRows(s_cond);
  auto result_condition_num = sigma.real().maxCoeff() / sigma.real().minCoeff();
  auto sigma_invsqrt = sigma.array().sqrt().inverse().matrix().asDiagonal();

  // make canonical X
  auto U_cond = U.block(0, s_rows - s_cond, s_rows, s_cond);
  X = U_cond * sigma_invsqrt;
  // make symmetric X
  if (symmetric) {
    X = X * U_cond.transpose().conjugate();
  }

  auto nbf_omitted = s_rows - s_cond;
  if (nbf_omitted < 0) throw "Error: dropping negative number of functions!";

  if (world.rank() == 0){
    std::cout << "\n\toverlap condition number = " << S_condition_num
              << " at k = " << k;
  }

  if (nbf_omitted > 0) {
      auto should_be_I = X.transpose().conjugate() * S_eig * X;
      auto I_real =
          Eigen::MatrixXd::Identity(should_be_I.rows(), should_be_I.cols());
      auto I_comp = I_real.template cast<std::complex<double>>();
      auto should_be_zero = (should_be_I - I_comp).norm();

      if (world.rank() == 0) {
        std::cout << " (dropped " << nbf_omitted << " "
                  << (nbf_omitted > 1 ? "fns" : "fn") << " to reduce to "
                  << result_condition_num << ")" << std::endl;
        std::cout << "\t\t||Xt*S*X - I||_2 = " << should_be_zero
                  << " (should be zero)" << std::endl;
      }

  }

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
std::vector<Matrixc> conditioned_orthogonalizer(
    const TArray overlap, int64_t k_size, double max_condition_num = 1.0e8) {
  std::vector<Matrixc> X;
  X.resize(k_size);

  auto tr0 = overlap.trange().data()[0];
  auto tr1 = overlap.trange().data()[1];

  assert(tr1.extent() == (tr0.extent() * k_size));

  auto overlap_eig = array_ops::array_to_eigen(overlap);
  for (auto k = 0; k < k_size; ++k) {
    std::vector<std::size_t> S_low{0, k * tr0.extent()};
    std::vector<std::size_t> S_up{tr0.extent(), (k + 1) * tr0.extent()};

    TArray S;
    S("mu, nu") = overlap("mu, nu").block(S_low, S_up);
    X[k] = gensqrtinv(S, false, max_condition_num, k);
  }

  return X;
}


}  // namespace scf
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_CONDITIONED_ORTHOGONALIZER_H_
