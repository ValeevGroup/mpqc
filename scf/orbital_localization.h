#pragma once
#ifndef MPQC_SCF_ORBITALLOCALIZATION_H
#define MPQC_SCF_ORBITALLOCALIZATION_H

#include "../include/tiledarray.h"
#include "../include/eigen.h"

#include "../ta_routines/array_to_eigen.h"
#include <cmath>


namespace mpqc {
namespace scf {

using Mat = Eig::Matrix<double, Eig::Dynamic, Eig::Dynamic, Eig::RowMajor>;

double norm(double x, double y, double z) {
    return std::sqrt(x * x + y * y + z * z);
}

Mat I_am_lazy_norm(Mat const &x, Mat const &y, Mat const &z) {
    Mat norms(x.rows(), x.cols());

    for (auto i = 0; i < x.rows(); ++i) {
        for (auto j = 0; j < x.cols(); ++j) {
            norms(i, j) = norm(x(i, j), y(i, j), z(i, j));
        }
    }
    return norms;
}

double boys_object(Mat const &x, Mat const &y, Mat const &z) {
    const auto x_sum = Eig::VectorXd(x.diagonal()).squaredNorm();
    const auto y_sum = Eig::VectorXd(y.diagonal()).squaredNorm();
    const auto z_sum = Eig::VectorXd(z.diagonal()).squaredNorm();
    return x_sum + y_sum + z_sum;
}

void calculate_iprs(Mat const &M, std::string title) {
    auto ipr = [](Eig::VectorXd const &v) {
        auto sum = 0.0;
        for (auto i = 0; i < v.size(); ++i) {
            sum += v[i] * v[i] * v[i] * v[i];
        }
        return std::pow(sum, 0.25);
    };

    auto largest_ipr = 0.0;
    auto smallest_ipr = 1.0;
    auto avg_ipr = 0.0;
    for (auto i = 0; i < M.cols(); ++i) {
        Eig::VectorXd col = M.col(i);
        col = 1 / col.norm() * col;
        auto val = ipr(col);
        avg_ipr += val;
        largest_ipr = std::max(largest_ipr, val);
        smallest_ipr = std::min(smallest_ipr, val);
    }
    std::cout << "Largest ipr in " << title << " = " << largest_ipr
              << std::endl;
    std::cout << "Smallest ipr in " << title << " = " << smallest_ipr
              << std::endl;
    std::cout << "Avg ipr in " << title << " = " << avg_ipr / M.cols()
              << std::endl;
}

Mat jacobi_sweeps(Mat const &C, std::vector<Mat> const &ao_xyz) {
    Mat mo_x = C.transpose() * ao_xyz[0] * C;
    Mat mo_y = C.transpose() * ao_xyz[1] * C;
    Mat mo_z = C.transpose() * ao_xyz[2] * C;

    auto gamma = [&](double Aij, double Bij) {
        auto ang = std::acos(-Aij / std::sqrt(Aij * Aij + Bij * Bij));
        if(ang > M_PI/2){
            ang = M_PI - ang;
        }  
        return ang/4;
    };

    Mat Cm = C;

    mo_x = Cm.transpose() * ao_xyz[0] * Cm;
    mo_y = Cm.transpose() * ao_xyz[1] * Cm;
    mo_z = Cm.transpose() * ao_xyz[2] * Cm;

    std::cout << "Input C = \n" << Cm << std::endl;
    auto iter = 1;
    while (iter < 3) {
        for (auto i = 0; i < Cm.cols(); ++i) {
            for (auto j = i + 1; j < Cm.cols(); ++j) {
                Vec3D vij = {mo_x(i,j), mo_y(i,j), mo_z(i,j)};
                Vec3D vi = {mo_x(i,i), mo_y(i,i), mo_z(i,i)};
                Vec3D vj = {mo_x(j,j), mo_y(j,j), mo_z(j,j)};
                Vec3D diff = vi - vj;

                auto nij = vij.norm();

                auto Aij = nij * nij - 0.25 * diff.squaredNorm();
                auto Bij = diff.dot(vij);// nij * diff.norm();

                auto g = gamma(Aij, Bij);
                auto cg = std::cos(g);
                auto sg = std::sin(g);

                Eig::VectorXd col_i = Cm.col(i);
                Eig::VectorXd col_j = Cm.col(j);

                Cm.col(i) = cg * col_i + sg * col_j;
                Cm.col(j) = cg * col_j - sg * col_i;

                if(Cm.cols() <= 5){
                    std::cout << "i j = " << i << " " << j << std::endl;
                    std::cout << "Cm = \n" << Cm << "\n" << std::endl;
                }

                mo_x = Cm.transpose() * ao_xyz[0] * Cm;
                mo_y = Cm.transpose() * ao_xyz[1] * Cm;
                mo_z = Cm.transpose() * ao_xyz[2] * Cm;
                auto crit = boys_object(mo_x, mo_y, mo_z);

                auto dm = I_am_lazy_norm(mo_x, mo_y, mo_z);
                Eig::SelfAdjointEigenSolver<Mat> es(dm);
                std::cout << "\tgamma = " << g * 4 << std::endl;
                std::cout << "\t\tCriteria = " << crit << std::endl;
                std::cout << "\t\tEvals dm = " << es.eigenvalues().transpose() << std::endl;
            }
        }
        std::cout << "Iteration " << iter++ << std::endl;
    }

    return Cm;
}


class BoysLocalization {
  public:
    template <typename Array>
    Array operator()(Array const &C, std::vector<Array> const &r_ao) const {

        auto c_eig = tcc::array_ops::array_to_eigen(C);

        auto ao_x = tcc::array_ops::array_to_eigen(r_ao[0]);
        auto ao_y = tcc::array_ops::array_to_eigen(r_ao[1]);
        auto ao_z = tcc::array_ops::array_to_eigen(r_ao[2]);

        auto cols = ao_x.cols();
        auto rows = ao_x.rows();

        Mat dm_eig;
        auto update_mo_dipole = [&]() {
            dm_eig = ao_x;
            for (auto i = 0; i < rows; ++i) {
                for (auto j = 0; j < cols; ++j) {
                    const auto x2 = ao_x(i, j) * ao_x(i, j);
                    const auto y2 = ao_y(i, j) * ao_y(i, j);
                    const auto z2 = ao_z(i, j) * ao_z(i, j);

                    dm_eig(i, j) = std::sqrt(x2 + y2 + z2);
                }
            }

            dm_eig = c_eig.transpose() * dm_eig * c_eig;
        };
        update_mo_dipole();

        // Eig::JacobiSVD<Mat> svd(dm_eig, Eig::ComputeThinU);
        // std::cout << "\nSingular values of dm = "
        //           << svd.singularValues().transpose() << std::endl;
        // std::cout << "Diagonal of dm " << dm_eig.diagonal().transpose()
        //           << std::endl;
        // std::cout << "sum squares evals = "
        //           << Eig::VectorXd(svd.singularValues()).squaredNorm()
        //           << std::endl;
        std::cout << "\n\nB func for input C "
                  << Eig::VectorXd(dm_eig.diagonal()).squaredNorm()
                  << std::endl;


        calculate_iprs(c_eig, "Input C");

        // Mat c_t = c_eig * svd.matrixU();
        // largest_ipr = 0.0;
        // smallest_ipr = 1.0;
        // avg_ipr = 0.0;
        // for (auto i = 0; i < c_t.cols(); ++i) {
        //     Eig::VectorXd col = c_t.col(i);
        //     col = 1 / col.norm() * col;
        //     auto val = ipr(col);
        //     avg_ipr += val;
        //     largest_ipr = std::max(largest_ipr, val);
        //     smallest_ipr = std::min(smallest_ipr, val);
        // }
        // std::cout << "\nLargest ipr in C trans = " << largest_ipr <<
        // std::endl;
        // std::cout << "Smallest ipr in C trans = " << smallest_ipr <<
        // std::endl;
        // std::cout << "Avg ipr in C trans = " << avg_ipr / c_t.cols()
        //           << std::endl;

        // if (c_t.cols() < 6) {
        //     std::cout << "C trans = \n" << c_t << std::endl;
        // }

        // auto update_mo_dipole_test = [&]() {
        //     dm_eig = ao_x;
        //     for (auto i = 0; i < rows; ++i) {
        //         for (auto j = 0; j < cols; ++j) {
        //             const auto x2 = ao_x(i, j) * ao_x(i, j);
        //             const auto y2 = ao_y(i, j) * ao_y(i, j);
        //             const auto z2 = ao_z(i, j) * ao_z(i, j);

        //             dm_eig(i, j) = std::sqrt(x2 + y2 + z2);
        //         }
        //     }

        //     dm_eig = c_t.transpose() * dm_eig * c_t;
        // };
        // update_mo_dipole_test();
        // std::cout << "sum square diag = "
        //           << Eig::VectorXd(dm_eig.diagonal()).squaredNorm()
        //           << std::endl;

        // test jacobi_sweeps
        std::cout << "\nStarting Boys" << std::endl;
        auto C_test = jacobi_sweeps(c_eig, {ao_x, ao_y, ao_z});
        calculate_iprs(C_test, "After Boys");

        // Doing silly thing
        Mat P = (c_eig.transpose() * c_eig).inverse() * c_eig.transpose();
        Mat U = P * Mat::Identity(c_eig.rows(), c_eig.cols());
        if(U.cols() < 6){
            std::cout << "\nInitial U = \n" << U << std::endl;
        }
        Eig::SelfAdjointEigenSolver<Mat> es(U * U.transpose());
        U = es.operatorInverseSqrt() * U;
        if(U.cols() < 6){
            std::cout << "Unitary U = \n" << U << "\n" << std::endl;
        }
        Eig::JacobiSVD<Mat> svdU(U);
        std::cout << "\nU Svals = " << svdU.singularValues().transpose() << std::endl;
        Mat C_exp = c_eig * U;
        calculate_iprs(C_exp, "Experimental C");
        if(C_exp.cols() < 6){
            std::cout << "My Expr C = \n" << C_exp << std::endl;
        }

        auto trange = C.trange();
        return tcc::array_ops::eigen_to_array<typename Array::value_type>(
              C.get_world(), C_test, trange.data()[0], trange.data()[1]);
    }
};

} // namespace scf
} // namespace mpqc

#endif //  MPQC_SCF_ORBITALLOCALIZATION_H
