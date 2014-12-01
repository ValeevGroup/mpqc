#include "../include/btas.h"
#include "../include/eigen.h"

#include <iostream>

int main(int argc, char *argv[]) {
    btas::Tensor<double> tensor(2, 4, 3);

    auto i = 1;
    tensor.generate([=]() mutable { return i++; });

    Eigen::MatrixXd mat1(2, 12);
    Eigen::MatrixXd mat2(4, 6);
    Eigen::MatrixXd mat3(3, 8);

    for (auto i = 0; i < 2; ++i) {
        std::cout << "Printing the " << i << " slice " << std::endl;
        for (auto j = 0; j < 4; ++j) {
            for (auto k = 0; k < 3; ++k) {
                std::cout << tensor(i, j, k) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }


    std::cout << "\nMatricization 1" << std::endl;
    for (auto i = 0; i < 2; ++i) {
        for (auto j = 0; j < 4; ++j) {
            for (auto k = 0; k < 3; ++k) {
                std::cout << tensor(i, j, k) << " ";
                mat1(i, j * 3 + k) = tensor(i, j, k);
            }
        }
        std::cout << std::endl;
    }

    std::cout << "\nMatricization 2" << std::endl;
    for (auto j = 0; j < 4; ++j) {
        for (auto i = 0; i < 2; ++i) {
            for (auto k = 0; k < 3; ++k) {
                std::cout << tensor(i, j, k) << " ";
                mat2(j, i * 3 + k) = tensor(i, j, k);
            }
        }
        std::cout << std::endl;
    }

    std::cout << "\nMatricization 3" << std::endl;
    for (auto k = 0; k < 3; ++k) {
        for (auto i = 0; i < 2; ++i) {
            for (auto j = 0; j < 4; ++j) {
                std::cout << tensor(i, j, k) << " ";
                mat3(k, i * 4 + j) = tensor(i, j, k);
            }
        }
        std::cout << std::endl;
    }

    std::cout << "\nEigen Mat 1 = \n" << mat1 << std::endl;
    std::cout << "Eigen Mat 2 = \n" << mat2 << std::endl;
    std::cout << "Eigen Mat 3 = \n" << mat3 << std::endl;


    Eigen::JacobiSVD<Eigen::MatrixXd> svd(mat1, Eigen::ComputeThinU
                                                | Eigen::ComputeThinV);
    Eigen::MatrixXd A1 = svd.matrixU().leftCols(2);
    std::cout << "SingularValues 1 = \n " << svd.singularValues().transpose()
              << std::endl;

    svd.compute(mat2, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd A2 = svd.matrixU().leftCols(2);
    std::cout << "SingularValues 2 = \n " << svd.singularValues().transpose()
              << std::endl;

    svd.compute(mat3, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd A3 = svd.matrixU().leftCols(2);
    std::cout << "SingularValues 3 = \n " << svd.singularValues().transpose()
              << std::endl;

    Eigen::MatrixXd Kron(A3.cols() * A2.cols(), A3.rows() * A2.rows());   
    std::cout << "Kron is " << Kron.rows() << " x " << Kron.cols() << std::endl;

    for(auto i = 0; i < A3.cols(); ++i){
        for(auto j = 0; j < A3.rows(); ++j){
            for(auto k = 0; k < A2.cols(); ++k){
                for(auto l = 0; l < A2.rows(); ++l){
                    Kron(i * A2.cols() + k, j * A2.rows() + l) = A3(j,i) * A2(l,k);
                }
            }
        }
    }

    std::cout << "Kron product = \n" << Kron << std::endl;
    std::cout << "Mat 1 dims = " << mat1.rows() << " x " << mat1.cols() << std::endl;

    Eigen::MatrixXd g_eig = A1.transpose() * mat1 * Kron.transpose();
    std::cout << "Finally G tensor = \n" << g_eig << std::endl;

    Eigen::MatrixXd mat1_approx = A1 * g_eig * Kron;
    std::cout << "Approximate mat1 = \n" << mat1_approx << std::endl;
    
    return 0;
}
