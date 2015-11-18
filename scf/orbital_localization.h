#pragma once
#ifndef MPQC_SCF_ORBITALLOCALIZATION_H
#define MPQC_SCF_ORBITALLOCALIZATION_H

#include "../include/tiledarray.h"
#include "../include/eigen.h"

#include "../ta_routines/array_to_eigen.h"
#include <cmath>

#include <array>
#include <iomanip>


namespace mpqc {
namespace scf {

using Mat = Eig::Matrix<double, Eig::Dynamic, Eig::Dynamic, Eig::RowMajor>;

double boys_object(std::array<Mat, 3> const &xyz) {
    auto sum = 0.0;
    for (auto i = 0; i < xyz[0].cols(); ++i) {
        sum += xyz[0](i, i) * xyz[0](i, i);
        sum += xyz[1](i, i) * xyz[1](i, i);
        sum += xyz[2](i, i) * xyz[2](i, i);
    }

    return sum;
    // const auto x_sum = Eig::VectorXd(xyz[0].diagonal()).squaredNorm();
    // const auto y_sum = Eig::VectorXd(xyz[1].diagonal()).squaredNorm();
    // const auto z_sum = Eig::VectorXd(xyz[2].diagonal()).squaredNorm();
    // return x_sum + y_sum + z_sum;
}

double gamma(double Aij, double Bij) {
    auto AB = std::sqrt(Aij * Aij + Bij * Bij);
    auto cos_gamma = -Aij / AB;
    auto sin_gamma = Bij / AB;
    auto ang = 0.25 * std::acos(cos_gamma) * ((sin_gamma < 0) ? -1 : 1);
    return (std::abs(ang) < 1e-7) ? 0 : ang;
};

void jacobi_sweeps(Mat &Cm, std::vector<Mat> const &ao_xyz) {
    std::array<Mat, 3> mo_xyz;
    mo_xyz[0] = Cm.transpose() * ao_xyz[0] * Cm;
    mo_xyz[1] = Cm.transpose() * ao_xyz[1] * Cm;
    mo_xyz[2] = Cm.transpose() * ao_xyz[2] * Cm;

    auto &mx = mo_xyz[0];
    auto &my = mo_xyz[1];
    auto &mz = mo_xyz[2];

    auto crit = boys_object(mo_xyz);
    auto iter = 1;
    auto error = crit - 0;
    while (error > 1e-6 && iter <= 150 ) {
        for (auto i = 0; i < Cm.cols(); ++i) {
            for (auto j = i + 1; j < Cm.cols(); ++j) {

                // double Aij = 0.0;
                // double Bij = 0.0;
                // for (auto z = 0; z < 3; ++z) {
                //     auto nij = mo_xyz[z](i, j);
                //     auto diff = mo_xyz[z](i, i) - mo_xyz[z](j, j);
                //     Aij += nij * nij - 0.25 * diff * diff;
                //     Bij += nij * diff;
                // }

                double Aij = mx(i, j) * mx(i, j) 
                           + my(i, j) * my(i, j)
                           + mz(i, j) * mz(i, j)
                           - 0.25 * (
                                   (mx(i, i) - mx(j, j)) 
                                 * (mx(i, i) - mx(j, j)) 
                                 + (my(i, i) - my(j, j))
                                 * (my(i, i) - my(j, j))
                                 + (mz(i, i) - mz(j, j))
                                 * (mz(i, i) - mz(j, j))
                              );
                double Bij = (mx(i,i) - mx(j,j)) * mx(i,j) 
                           + (my(i,i) - my(j,j)) * my(i,j) 
                           + (mz(i,i) - mz(j,j)) * mz(i,j);

                auto g = gamma(Aij, Bij);
                auto cg = std::cos(g);
                auto sg = std::sin(g);

                Eig::VectorXd col_i = Cm.col(i);
                Eig::VectorXd col_j = Cm.col(j);

                Cm.col(i) = cg * col_i + sg * col_j;
                Cm.col(j) = -sg * col_i + cg * col_j;

                for (auto z = 0; z < 3; ++z) {
                    auto &m = mo_xyz[z];
                    Eig::VectorXd z_i = m.col(i);
                    Eig::VectorXd z_j = m.col(j);

                    m.col(i) = cg * z_i + sg * z_j;
                    m.col(j) = -sg * z_i + cg * z_j;

                    z_i = m.row(i);
                    z_j = m.row(j);

                    m.row(i) = cg * z_i + sg * z_j;
                    m.row(j) = -sg * z_i + cg * z_j;
                }
            }
        }
        auto old_crit = crit;
        crit = boys_object(mo_xyz);
        error = std::abs(old_crit - crit);
        std::cout << "Iteration: " << iter++ << " cost: " << crit << " error "
                  << error << std::endl;
    }
}


class BoysLocalization {
  public:
    template <typename Array>
    Array operator()(Array const &C, std::vector<Array> const &r_ao) const {
        auto ao_x = tcc::array_ops::array_to_eigen(r_ao[0]);
        auto ao_y = tcc::array_ops::array_to_eigen(r_ao[1]);
        auto ao_z = tcc::array_ops::array_to_eigen(r_ao[2]);
        auto c_eig = tcc::array_ops::array_to_eigen(C);

        Mat mo_x = c_eig.transpose() * ao_x * c_eig;
        Mat mo_y = c_eig.transpose() * ao_y * c_eig;
        Mat mo_z = c_eig.transpose() * ao_z * c_eig;


        std::cout << "\nStarting Boys" << std::endl;
        jacobi_sweeps(c_eig, {ao_x, ao_y, ao_z});

        mo_x = c_eig.transpose() * ao_x * c_eig;
        mo_y = c_eig.transpose() * ao_y * c_eig;
        mo_z = c_eig.transpose() * ao_z * c_eig;

        auto trange = C.trange();
        return tcc::array_ops::eigen_to_array<typename Array::value_type>(
              C.get_world(), c_eig, trange.data()[0], trange.data()[1]);
    }
};

class StuffFromFactorAna {
  public:
    template <typename Array>
    Array operator()(Array const &C) const {
        auto A = tcc::array_ops::array_to_eigen(C);

        auto errorQ = [](Mat const &L) {
            auto sum = 0.0;
            for (auto i = 0; i < L.rows(); ++i) {
                for (auto j = 0; j < L.cols(); ++j) {
                    auto val = L(i, j);
                    sum += val * val * val * val;
                }
            }

            return 0.25 * sum;
        };

        if (A.cols() < 5) {
            std::cout << std::setprecision(3);
            std::cout << "Input A = \n" << A << std::endl;
            std::cout << std::setprecision(15);
        }
        Mat T = Mat::Identity(A.cols(), A.cols());
        Mat L = A * T;
        auto max = errorQ(L);
        auto max_old = 0.0;
        auto error = max - max_old;
        auto iter = 0;
        while (iter < 30 && error > 1e-2) {
            ++iter;
            L = L.array() * L.array() * L.array();
            Mat G = A.transpose() * (L);
            Eig::JacobiSVD<Mat> svd(G + 0.5 * T,
                                    Eig::ComputeThinU | Eig::ComputeThinV);
            T = svd.matrixU() * svd.matrixV().transpose();
            L = A * T;
            max_old = max;
            max = errorQ(L);
            error = max - max_old;
            std::cout << "Iter: " << iter << " Error in Q = " << error
                      << std::endl;
        }
        A = A * T;
        if (A.cols() < 5) {
            std::cout << std::setprecision(3);
            std::cout << "Output A = \n" << A << std::endl;
            std::cout << std::setprecision(15);
        }

        auto trange = C.trange();
        return tcc::array_ops::eigen_to_array<typename Array::value_type>(
              C.get_world(), A, trange.data()[0], trange.data()[1]);
    }
};

} // namespace scf
} // namespace mpqc

#endif //  MPQC_SCF_ORBITALLOCALIZATION_H
