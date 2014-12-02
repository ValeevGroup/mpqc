#include "../include/btas.h"
#include "../include/eigen.h"

#include <iostream>

Eigen::MatrixXd matricization(btas::Tensor<double> const &input, int i) {
    if (i == 0) {
        for(auto y : input.extent()){
            std::cout << y << std::endl;
        }
        for(auto it = input.extent().begin(); it != input.extent().end(); ++it){
            std::cout << *it << std::endl;
        }

        auto cols = 1ul;
        for(auto x : input.extent()){
            cols *= x;
        }
        cols /= input.extent().front();


        Eigen::MatrixXd mat(input.extent().front(), cols);
        std::cout << "Mat dims = " << mat.rows() << " x " << mat.cols()
                  << std::endl;

        std::memcpy(mat.data(), input.storage().data(), mat.size());

        return mat;
    }

    return Eigen::MatrixXd::Zero(i,i);
}

void tucker_decomp(btas::Tensor<double> const &input) {
    std::vector<Eigen::MatrixXd> mats(input.rank());

    mats[0] = matricization(input, 0);
}

int main() {
    btas::Tensor<double> tensor(2, 4, 3);

    auto i = 1;
    tensor.generate([=]() mutable { return i++; });

    Eigen::MatrixXd mat1(2, 12);
    Eigen::MatrixXd mat2(4, 6);
    Eigen::MatrixXd mat3(3, 8);


    for (auto i = 0; i < 2; ++i) {
        for (auto j = 0; j < 4; ++j) {
            for (auto k = 0; k < 3; ++k) {
                mat1(i, j * 3 + k) = tensor(i, j, k);
            }
        }
    }

    for (auto j = 0; j < 4; ++j) {
        for (auto i = 0; i < 2; ++i) {
            for (auto k = 0; k < 3; ++k) {
                mat2(j, i * 3 + k) = tensor(i, j, k);
            }
        }
    }

    for (auto k = 0; k < 3; ++k) {
        for (auto i = 0; i < 2; ++i) {
            for (auto j = 0; j < 4; ++j) {
                mat3(k, i * 4 + j) = tensor(i, j, k);
            }
        }
    }

    std::cout << "\nEigen Mat 1 = \n" << mat1 << std::endl;
    std::cout << "Eigen Mat 2 = \n" << mat2 << std::endl;
    std::cout << "Eigen Mat 3 = \n" << mat3 << std::endl;


    Eigen::JacobiSVD
        <Eigen::MatrixXd> svd(mat1, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd A1 = svd.matrixU().leftCols(2);

    svd.compute(mat2, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd A2 = svd.matrixU().leftCols(2);

    svd.compute(mat3, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd A3 = svd.matrixU().leftCols(2);

    Eigen::MatrixXd Kron(A3.cols() * A2.cols(), A3.rows() * A2.rows());

    for (auto i = 0; i < A3.cols(); ++i) {
        for (auto j = 0; j < A3.rows(); ++j) {
            for (auto k = 0; k < A2.cols(); ++k) {
                for (auto l = 0; l < A2.rows(); ++l) {
                    Kron(i * A2.cols() + k, j *A2.rows() + l) = A3(j, i)
                                                                * A2(l, k);
                }
            }
        }
    }


    Eigen::MatrixXd g_eig = A1.transpose() * mat1 * Kron.transpose();

    Eigen::MatrixXd mat1_approx = A1 * g_eig * Kron;
    std::cout << "Approximate mat1 = \n" << mat1_approx << std::endl;

    tucker_decomp(tensor);

    return 0;
}
