
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

class BoysLocalization {
 public:
  template <typename Array>
  Array operator()(Array const &C, std::vector<Array> const &r_ao) const {
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
