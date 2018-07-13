
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_ORBITAL_LOCALIZATION_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_ORBITAL_LOCALIZATION_H_

#include <cmath>
#include <array>
#include <iomanip>

#include <tiledarray.h>

#include "mpqc/math/external/eigen/eigen.h"
#include "mpqc/math/tensor/clr/array_to_eigen.h"
#include "mpqc/util/keyval/keyval.h"
#include "mpqc/util/misc/assert.h"
#include "mpqc/util/core/exception.h"

namespace mpqc {
namespace lcao {

/// Localizes orbitals using LCAO-specific info (e.g. AO-basis operators)
template <typename Tile, typename Policy>
class OrbitalLocalizer : public DescribedClass {
 public:
  /// @param[in] C input LCAOs
  /// @param[in] ncols_of_C_to_skip the number of columns of C to keep
  /// non-localized
  /// @return transformation matrix U that converts C to localized LCAOs
  virtual TA::DistArray<Tile, Policy> compute(
      TA::DistArray<Tile, Policy> const &C,
      size_t ncols_of_C_to_skip = 0) const = 0;

  /// this or its Eigen counterpart must be called before compute()
  OrbitalLocalizer &initialize(TA::DistArray<Tile, Policy> S_ao,
                               std::vector<TA::DistArray<Tile, Policy>> mu_ao) {
    ao_s_ = math::array_to_eigen(S_ao);
    ao_x_ = math::array_to_eigen(mu_ao[0]);
    ao_y_ = math::array_to_eigen(mu_ao[1]);
    ao_z_ = math::array_to_eigen(mu_ao[2]);
    initialized_ = true;
    return *this;
  }

  /// this or its TA counterpart must be called before compute()
  OrbitalLocalizer &initialize(const math::Matrix<typename Tile::value_type>& S_ao,
                               const std::vector<math::Matrix<typename Tile::value_type>>& mu_ao) {
    ao_s_ = S_ao;
    ao_x_ = mu_ao[0];
    ao_y_ = mu_ao[1];
    ao_z_ = mu_ao[2];
    initialized_ = true;
    return *this;
  }

 protected:
  math::Matrix<typename Tile::value_type> ao_s_, ao_x_, ao_y_, ao_z_;
  bool initialized_ = false;
};

namespace scf {

using Mat =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

/// finds a stationary point of the Foster-Boys objective function
/// @param[in,out] Cm on input: LCAOs to be localized; on output: localized
/// LCAOs
/// @param[in,out] U on output: transformation matrix converting original to
/// localized LCAOs
/// @param[in] ao_xyz the {x,y,z} electric dipole integral matrices in AO basis
/// @param convergence_threshold stop once maximum rotation angle (in rad)
/// changes between iterations by less than this
/// @param max_iter do not exceed this many iterations
/// @return true if converged
bool fb_jacobi_sweeps(Mat &Cm, Mat &U, std::vector<Mat> const &ao_xyz,
                      double convergence_threshold, size_t max_iter);

/// Performs Foster-Boys localization
/// (see J. Foster and S. Boys, Rev Mod Phys 32, 300 (1960)).
template <typename Tile, typename Policy>
class FosterBoysLocalizer : public OrbitalLocalizer<Tile, Policy> {
 public:
  // clang-format off
  /**
   * KeyVal constructor for FosterBoysLocalizer
   *
   * @param kv the KeyVal object; it will be queried for all keywords of AOWavefunction as well as the following additional keywords:
   * | Keyword | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | @c convergence | double | 1e-8 | the Jacobi solver is converged when the maximum rotation angle (in rad) does not exceed this value  |
   * | @c max_iter | int | 50 | the maximum number of Jacobi iterations |
   */
  // clang-format on
  FosterBoysLocalizer(const KeyVal &kv = KeyVal{})
      : jacobi_convergence_threshold_(kv.value<double>("convergence", 1e-8, KeyVal::is_nonnegative)),
        jacobi_max_iter_(kv.value<std::size_t>("max_iter", 50)) {}

  /// @param C input LCAOs
  /// @param {x,y,z} electric dipole operator matrices, in AO basis
  /// @param[in] ncols_of_C_to_skip the number of columns of C to keep
  ///            non-localized, presumably because they are already localized
  /// @return transformation matrix U that converts C to localized LCAOs
  TA::DistArray<Tile, Policy> compute(
      TA::DistArray<Tile, Policy> const &C,
      size_t ncols_of_C_to_skip = 0) const override {
    MPQC_ASSERT(this->initialized_);
    auto c_eig = math::array_to_eigen(C);

    auto trange = C.trange();
    return math::eigen_to_array<Tile, Policy>(
        C.world(), (*this)(c_eig, this->ao_x_, this->ao_y_, this->ao_z_, ncols_of_C_to_skip),
        trange.data()[1], trange.data()[1]);
  }

 private:
  /// @param[in,out] C on input: LCAO coefficients, on output: localized LCAO
  /// coefficients
  /// @param {x,y,z} electric dipole operator matrices, in AO basis
  /// @param[in] ncols_of_C_to_skip the number of columns of C to keep
  ///            non-localized, presumably because they are already localized
  /// @return transformation matrix U that converts C to localized LCAOs, i.e.
  /// \c Cao("mu,k") * U("k,i") computes the AO coefficients of localized MOs
  /// from the AO coefficients of input MOs";
  template <typename EigMat>
  EigMat operator()(EigMat &C, const EigMat &ao_x, const EigMat &ao_y,
                    const EigMat &ao_z, size_t ncols_of_C_to_skip = 0) const {
    EigMat C_loc =
        C.block(0, ncols_of_C_to_skip, C.rows(), C.cols() - ncols_of_C_to_skip);

    EigMat U_loc = EigMat::Identity(C_loc.cols(), C_loc.cols());
    auto converged = fb_jacobi_sweeps(C_loc, U_loc, {ao_x, ao_y, ao_z},
                                      jacobi_convergence_threshold_, jacobi_max_iter_);
    if (!converged) {
      std::ostringstream oss;
      oss << "Foster-Boys Jacobi sweeps failed to converge to threshold=" << jacobi_convergence_threshold_
          << " in " << jacobi_max_iter_ << " iterations";
      throw AlgorithmException(oss.str().c_str(), __FILE__, __LINE__);
    }

    EigMat U = EigMat::Identity(C.cols(), C.cols());
    U.block(ncols_of_C_to_skip, ncols_of_C_to_skip, U_loc.rows(),
            U_loc.cols()) = U_loc;
    return U;
  }

  double jacobi_convergence_threshold_;
  size_t jacobi_max_iter_;
};

/// Performs Rank Revealing QR localization
/// (see Damle, A. et al., J. Chem. Theory Comput. 2015, 11 (4), 1463.)
template <typename Tile, typename Policy>
class RRQRLocalizer : public OrbitalLocalizer<Tile,Policy> {
 public:
  RRQRLocalizer(const KeyVal &) {}

  TA::DistArray<Tile, Policy> compute(TA::DistArray<Tile, Policy> const &C,
                                      size_t ncols_of_C_to_skip = 0) const override {
    MPQC_ASSERT(this->initialized_);
    auto c_eig = math::array_to_eigen(C);
    auto trange = C.trange();
    return math::eigen_to_array<Tile, Policy>(
        C.world(), (*this)(c_eig, this->ao_s_, ncols_of_C_to_skip), trange.data()[1],
        trange.data()[1]);
  }

  /// @param[in, out] C on input: LCAO coefficients, on output: localizeed LCAO
  /// coefficients
  /// @param[in] S Overlap matrix
  /// @param[in] ncols_of_C_to_skip the number of columns of C to keep
  /// non-localized
  /// @return transformation matrix U that converts C to localized LCAOs, i.e.
  /// \c Cao("mu,k") * U("k,i") computes the AO coefficients of localized MOs
  /// from the AO coefficients of input MOs";
  template <typename EigMat>
  EigMat operator()(EigMat &C, const EigMat &S,
                    size_t ncols_of_C_to_skip = 0) const {
    EigMat C_loc =
        C.block(0, ncols_of_C_to_skip, C.rows(), C.cols() - ncols_of_C_to_skip);

    Eigen::SelfAdjointEigenSolver<EigMat> eig_solver(S);
    EigMat X = eig_solver.operatorInverseSqrt();
    EigMat Psi_tr = ((X.inverse()) * C_loc).transpose();

    Eigen::ColPivHouseholderQR<RowMatrixXd> qr(Psi_tr);
    EigMat Q = qr.householderQ().setLength(qr.nonzeroPivots());

    EigMat U = Eigen::MatrixXd::Identity(C.cols(), C.cols());
    U.block(ncols_of_C_to_skip, ncols_of_C_to_skip, Q.rows(), Q.cols()) = Q;
    return U;
  }
};

}  // namespace scf
}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_ORBITAL_LOCALIZATION_H_
