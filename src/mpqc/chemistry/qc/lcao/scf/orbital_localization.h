
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

double boys_object(std::array<Mat, 3> const &xyz);

double gamma(double Aij, double Bij);
void jacobi_sweeps(Mat &Cm, Mat &U, std::vector<Mat> const &ao_xyz);

/// Performs Boys-Foster localization
class BoysLocalization {
 public:
  /// @param[in] C input vectors (rows = AOs)
  /// @param[in] r_ao input 1st (dipole) moment vectors, in AO basis
  /// @param[in] ncols_of_C_to_skip the number of columns of C to keep
  ///            non-localized, presumably because they are already localized
  /// @return unitary transformation that "localizes" C
  template <typename Tile, typename Policy>
  TA::DistArray<Tile,Policy> operator()(TA::DistArray<Tile,Policy> const &C,
                                        std::vector<TA::DistArray<Tile,Policy>> const &r_ao,
                                        size_t ncols_of_C_to_skip = 0) const {
    auto ao_x = array_ops::array_to_eigen(r_ao[0]);
    auto ao_y = array_ops::array_to_eigen(r_ao[1]);
    auto ao_z = array_ops::array_to_eigen(r_ao[2]);
    auto c_eig = array_ops::array_to_eigen(C);
    decltype(c_eig) c_eig_loc = c_eig.block(0 , ncols_of_C_to_skip, c_eig.rows(), c_eig.cols() - ncols_of_C_to_skip);

    RowMatrixXd U_loc = RowMatrixXd::Identity(c_eig_loc.cols(), c_eig_loc.cols());
    jacobi_sweeps(c_eig_loc, U_loc, {ao_x, ao_y, ao_z});

    RowMatrixXd U = RowMatrixXd::Identity(c_eig.cols(), c_eig.cols());
    U.block(ncols_of_C_to_skip, ncols_of_C_to_skip, U_loc.rows(), U_loc.cols()) = U_loc;

    auto trange = C.trange();
    return array_ops::eigen_to_array<Tile,Policy>(
        C.world(), U, trange.data()[1], trange.data()[1]);
  }
};

}  // namespace scf
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_ORBITAL_LOCALIZATION_H_
