
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
  /// @return transformation matrix U that converts C to localized LCAOs
  template <typename Tile, typename Policy>
  TA::DistArray<Tile,Policy> operator()(TA::DistArray<Tile,Policy> const &C, std::vector<TA::DistArray<Tile,Policy>> const &r_ao) const {
    auto ao_x = array_ops::array_to_eigen(r_ao[0]);
    auto ao_y = array_ops::array_to_eigen(r_ao[1]);
    auto ao_z = array_ops::array_to_eigen(r_ao[2]);
    auto c_eig = array_ops::array_to_eigen(C);

    RowMatrixXd U = RowMatrixXd::Identity(c_eig.cols(), c_eig.cols());
    jacobi_sweeps(c_eig, U, {ao_x, ao_y, ao_z});

    auto trange = C.trange();
    return array_ops::eigen_to_array<Tile,Policy>(
        C.world(), U, trange.data()[1], trange.data()[1]);
  }
};

}  // namespace scf
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_ORBITAL_LOCALIZATION_H_
