
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_ORBITAL_LOCALIZATION_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_ORBITAL_LOCALIZATION_H_

#include "mpqc/math/external/eigen/eigen.h"
#include <tiledarray.h>

#include <cmath>

#include <array>
#include <iomanip>

#include "mpqc/math/external/eigen/eigen.h"
#include "mpqc/math/tensor/clr/array_to_eigen.h"

namespace mpqc {
namespace scf {

using Mat =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

/// @param[in,out] Cm on input: LCAOs to be localized; on output: localized LCAOs
/// @param[in,out] U on output: transformation matrix converting original to localized LCAOs
/// @param[in] ao_xyz the {x,y,z} electric dipole integral matrices in AO basis
/// @param convergence_threshold stop once maximum rotation angle (in rad) changes between iterations by less than this
/// @param max_iter do not exceed this many iterations
void jacobi_sweeps(Mat &Cm, Mat &U, std::vector<Mat> const &ao_xyz,
                   double convergence_threshold = 1e-4,
                   size_t max_iter = 50);

/// Performs Foster-Boys localization (see J. Foster and S. Boys, Rev Mod Phys 32, 300 (1960)).
class FosterBoysLocalization {
 public:
  /// @param C input LCAOs
  /// @param {x,y,z} electric dipole operator matrices, in AO basis
  /// @param[in] ncols_of_C_to_skip the number of columns of C to keep
  ///            non-localized, presumably because they are already localized
  /// @return transformation matrix U that converts C to localized LCAOs
  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> operator()(
      TA::DistArray<Tile, Policy> const &C,
      std::vector<TA::DistArray<Tile, Policy>> const &r_ao,
      size_t ncols_of_C_to_skip = 0) const {
    auto ao_x = array_ops::array_to_eigen(r_ao[0]);
    auto ao_y = array_ops::array_to_eigen(r_ao[1]);
    auto ao_z = array_ops::array_to_eigen(r_ao[2]);
    auto c_eig = array_ops::array_to_eigen(C);

    auto trange = C.trange();
    return array_ops::eigen_to_array<Tile, Policy>(
        C.world(), (*this)(c_eig, ao_x, ao_y, ao_z, ncols_of_C_to_skip),
        trange.data()[1], trange.data()[1]);
  }

  /// @param[in,out] C on input: LCAO coeffcients, on output: localized LCAO coefficients
  /// @param {x,y,z} electric dipole operator matrices, in AO basis
  /// @param[in] ncols_of_C_to_skip the number of columns of C to keep
  ///            non-localized, presumably because they are already localized
  /// @return transformation matrix U that converts C to localized LCAOs
  /// U is a matrix such that: Cao("mu,i") = Cao("mu,k") * U("k,i");
  template <typename EigMat>
  EigMat operator()(EigMat &C, const EigMat &ao_x, const EigMat &ao_y,
                    const EigMat &ao_z, size_t ncols_of_C_to_skip = 0) const {
    EigMat C_loc =
        C.block(0, ncols_of_C_to_skip, C.rows(), C.cols() - ncols_of_C_to_skip);

    EigMat U_loc = EigMat::Identity(C_loc.cols(), C_loc.cols());
    jacobi_sweeps(C_loc, U_loc, {ao_x, ao_y, ao_z});

    EigMat U = EigMat::Identity(C.cols(), C.cols());
    U.block(ncols_of_C_to_skip, ncols_of_C_to_skip, U_loc.rows(),
            U_loc.cols()) = U_loc;
    return U;
  }
};

/// Performs Rank Revealing QR localization
class RRQRLocalization{
public:
  /// @param[in] C on input: LCAO coefficients
  /// @param[in] S : Overlap matrix
  /// @param[in] ncols_of_C_to_skip the number of columns of C to keep non-localized
  /// @return  transformation matrix U that converts C to localized LCAOs
  /// U is a matrix such that: Cao("mu,i") = Cao("mu,k") * U("k,i");
    template <typename Tile, typename Policy>
    TA::DistArray<Tile, Policy> operator()(
            TA::DistArray<Tile, Policy> &C,
            TA::DistArray<Tile, Policy> const &S,
            size_t ncols_of_C_to_skip = 0) const {

        auto c_eig = array_ops::array_to_eigen(C);
        auto s_eig = array_ops::array_to_eigen(S);
        auto trange = C.trange();
        return array_ops::eigen_to_array<Tile, Policy>(
                C.world(), (*this)(c_eig, s_eig, ncols_of_C_to_skip),
	            trange.data()[1], trange.data()[1]);
    }

    template <typename EigMat>
	EigMat operator()(EigMat& C, const EigMat& S,
                  size_t ncols_of_C_to_skip = 0) const {
    EigMat C_loc =
            C.block(0, ncols_of_C_to_skip, C.rows(), C.cols() - ncols_of_C_to_skip);

    Eigen::SelfAdjointEigenSolver<EigMat> eig_solver(S);
    EigMat X = eig_solver.operatorInverseSqrt();
    EigMat Psi_tr = ((X.inverse())*C_loc).transpose();

    Eigen::ColPivHouseholderQR<RowMatrixXd> qr(Psi_tr);
    EigMat Q = qr.householderQ().setLength(qr.nonzeroPivots());

    EigMat U = Eigen::MatrixXd::Identity(C.cols(),C.cols());
    U.block(ncols_of_C_to_skip, ncols_of_C_to_skip, Q.rows(), Q.cols()) = Q;
	return U;
  }
};

}  // namespace scf
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_ORBITAL_LOCALIZATION_H_
