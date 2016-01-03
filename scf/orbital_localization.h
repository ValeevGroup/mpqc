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
}

double gamma(double Aij, double Bij) {
    auto AB = std::sqrt(Aij * Aij + Bij * Bij);
    auto cos_gamma = -Aij / AB;
    auto sin_gamma = Bij / AB;
    auto ang = 0.25 * std::acos(cos_gamma) * ((sin_gamma < 0) ? -1 : 1);
    return (std::abs(ang) < 1e-7) ? 0 : ang;
};

void jacobi_sweeps(Mat &Cm, Mat &U, std::vector<Mat> const &ao_xyz) {
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
    while (error > 1e-8 && iter <= 150) {
        for (auto i = 0; i < Cm.cols(); ++i) {
            for (auto j = i + 1; j < Cm.cols(); ++j) {

                Vec3D vij = {mx(i, j), my(i, j), mz(i, j)};
                Vec3D vii = {mx(i, i), my(i, i), mz(i, i)};
                Vec3D vjj = {mx(j, j), my(j, j), mz(j, j)};

                double Aij = vij.squaredNorm()
                             - 0.25 * (vii - vjj).squaredNorm();
                double Bij = (vii - vjj).dot(vij);

                auto g = gamma(Aij, Bij);
                auto cg = std::cos(g);
                auto sg = std::sin(g);

                Eig::VectorXd col_Ui = U.col(i);
                Eig::VectorXd col_Uj = U.col(j);

                U.col(i) = cg * col_Ui + sg * col_Uj;
                U.col(j) = -sg * col_Ui + cg * col_Uj;

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
        ++iter;
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

        // Try an approximate initial guess
        MatrixD U = MatrixD::Identity(c_eig.cols(), c_eig.cols());

        jacobi_sweeps(c_eig, U, {ao_x, ao_y, ao_z});

        auto trange = C.trange();
        return tcc::array_ops::eigen_to_array<typename Array::value_type>(
              C.get_world(), U, trange.data()[1], trange.data()[1]);
    }
};

// Taken from APILPraat/dwtools/Configuration.cpp on Github under GPL2
class JacobiVarimax {
  public:
    double varimax_object(Mat const &m) {
        auto v4 = 0.0;
        for (auto j = 0; j < m.cols(); ++j) {
            double sum4 = 0.0, mean = 0.0;

            for (auto i = 0; i < m.rows(); ++i) {
                double sq = m(i, j) * m(i, j);
                sum4 += sq * sq;
                mean += sq;
            }
            v4 += sum4;
            v4 -= (mean * mean) / m.rows();
        }

        return v4;
    }

    template <typename Array>
    Array operator()(Array const &C) {
        Mat A = tcc::array_ops::array_to_eigen(C);

        auto object = varimax_object(A);
        std::cout << "Initial objective of C = " << object << std::endl;
        auto error = object;
        auto max_iter = 100;
        auto iter = 0;
        while (iter < max_iter && error > 1e-8) {
            for (auto x = 0; x < A.cols(); ++x) {
                for (auto y = 0; y < x; ++y) {
                    Eig::VectorXd vx = A.col(x).array() * A.col(x).array();
                    Eig::VectorXd vy = A.col(y).array() * A.col(y).array();

                    Eig::VectorXd u = vx - vy;
                    Eig::VectorXd v = 2 * vx.array() * vy.array();

                    auto uavg = u.array().sum() / u.size();
                    auto vavg = v.array().sum() / v.size();

                    u = u.array() - uavg;
                    v = u.array() - vavg;

                    double a = (u.array() * v.array()).sum();
                    Eig::VectorXd u2 = u.array() * u.array();
                    Eig::VectorXd v2 = u.array() * u.array();
                    double b = (u2 - v2).sum();

                    double c = std::sqrt(4 * a * a + b * b);
                    double w = std::sqrt((c + b) / (2 * c));
                    if (a > 0) {
                        w = -w;
                    }

                    double cost = std::sqrt(0.5 + 0.5 * w);
                    double sint = std::sqrt(0.5 - 0.5 * w);
                    if (std::acos(cost) != std::asin(sint)) {
                        std::cout << "Phi disc = " << std::acos(cost) << " "
                                  << std::asin(sint) << std::endl;
                    }
                    if (std::acos(cost) < 1e-7) {
                        continue;
                    }

                    double t22 = cost;
                    double t11 = cost;
                    double t12 = -sint;
                    double t21 = sint;
                    if (w < 0) {
                        t11 = sint;
                        t12 = t21 = cost;
                        t22 = -sint;
                    }

                    Eig::VectorXd col_x = A.col(x);
                    Eig::VectorXd col_y = A.col(y);

                    auto a_norm = A.lpNorm<2>();
                    A.col(x) = t11 * col_x + t21 * col_y;
                    A.col(y) = t12 * col_x + t22 * col_y;
                    auto a_new_norm = A.lpNorm<2>();
                    if (std::abs(a_norm - a_new_norm) > 1e-15) {
                        std::cout << "Norm fucked up by rotation!" << std::endl;
                    }
                }
            }
            auto old_object = object;
            object = varimax_object(A);
            error = std::abs(object - old_object);
            std::cout << "Iter: " << ++iter << " varimax object: " << object
                      << " error: " << error << std::endl;
        }

        auto trange = C.trange();
        return tcc::array_ops::eigen_to_array<typename Array::value_type>(
              C.get_world(), A, trange.data()[0], trange.data()[1]);
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

        Mat T = Mat::Identity(A.cols(), A.cols());
        Mat L = A * T;
        auto max = errorQ(L);
        auto max_old = 0.0;
        auto error = max - max_old;
        auto iter = 0;
        double alpha = 0;
        double dalpha = 0.5;
        while (iter < 40 && error > 1e-6) {
            ++iter;
            L = L.array() * L.array() * L.array();
            Mat G = A.transpose() * (L);
            if (max > 0.1) {
                Eig::JacobiSVD<Mat> svd(G,
                                        Eig::ComputeThinU | Eig::ComputeThinV);
                T = svd.matrixU() * svd.matrixV().transpose();
            } else {
                alpha += dalpha;
                Eig::JacobiSVD<Mat> svd(G + alpha * T,
                                        Eig::ComputeThinU | Eig::ComputeThinV);
                T = svd.matrixU() * svd.matrixV().transpose();
            }
            L = A * T;
            max_old = max;
            max = errorQ(L);
            error = max - max_old;
            std::cout << "Iter: " << iter << " Error in Q = " << error
                      << std::endl;
        }
        A = A * T;

        auto trange = C.trange();
        return tcc::array_ops::eigen_to_array<typename Array::value_type>(
              C.get_world(), A, trange.data()[0], trange.data()[1]);
    }
};

} // namespace scf
} // namespace mpqc

#endif //  MPQC_SCF_ORBITALLOCALIZATION_H
