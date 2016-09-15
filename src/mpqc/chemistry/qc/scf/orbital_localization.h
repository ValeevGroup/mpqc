#pragma once
#ifndef MPQC_SCF_ORBITALLOCALIZATION_H
#define MPQC_SCF_ORBITALLOCALIZATION_H

#include "../../../../../include/tiledarray.h"
#include "../../../../../include/eigen.h"

#include "../../../../../ta_routines/array_to_eigen.h"
#include <cmath>

#include <array>
#include <iomanip>


namespace mpqc {
namespace scf {

using Mat = Eig::Matrix<double, Eig::Dynamic, Eig::Dynamic, Eig::RowMajor>;

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

        MatrixD U = MatrixD::Identity(c_eig.cols(), c_eig.cols());
        jacobi_sweeps(c_eig, U, {ao_x, ao_y, ao_z});

        auto trange = C.trange();
        return array_ops::eigen_to_array<typename Array::value_type>(
              C.get_world(), U, trange.data()[1], trange.data()[1]);
    }
};

} // namespace scf
} // namespace mpqc

#endif //  MPQC_SCF_ORBITALLOCALIZATION_H
