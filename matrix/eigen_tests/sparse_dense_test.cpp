#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <random>
#include <chrono>
#include <string>
#include <madness/tensor/cblas.h>

using namespace Eigen;

int main(int argc, char **argv) {
    int input = (argc > 1) ? std::stoi(argv[1]) : 1000;
    int num_insig = (argc > 2) ? std::stoi(argv[2]) : 990;

    MatrixXd mat = MatrixXd::Random(input, input);
    std::mt19937 gen(100);
    std::uniform_int_distribution<int> dist(0, input - 1);

    for (auto j = 0; j < mat.cols(); ++j) {
        for (auto i = 0; i < num_insig; ++i) {
            mat(j, dist(gen)) = 0.0;
        }
    }

    if (mat.size() < 100) {
        std::cout << "mat = \n" << mat << std::endl;
    }

    long num_zeros = 0;
    for (auto i = 0; i < mat.cols(); ++i) {
        for (auto j = 0; j < mat.rows(); ++j) {
            num_zeros += (mat(i, j) == 0) ? 1 : 0;
        }
    }

    MatrixXd A = mat;
    MatrixXd B = mat;
    MatrixXd C = mat;

    double percent_zero = double(num_zeros) / double(mat.size()) * 100;


    SparseMatrix<double> sA = mat.sparseView();
    SparseMatrix<double> sB = mat.sparseView();
    SparseMatrix<double> sC = mat.sparseView();
    sA.makeCompressed();
    sB.makeCompressed();
    sC.makeCompressed();

    auto sparse_mult0 = std::chrono::steady_clock::now();
    sC = sA * sB + sC;
    auto sparse_mult1 = std::chrono::steady_clock::now();
    double sparse_time
        = std::chrono::duration_cast
          <std::chrono::duration<double>>(sparse_mult1 - sparse_mult0).count();

    auto dense_mult0 = std::chrono::steady_clock::now();
    C = A * B + C;
    auto dense_mult1 = std::chrono::steady_clock::now();
    double dense_time
        = std::chrono::duration_cast
          <std::chrono::duration<double>>(dense_mult1 - dense_mult0).count();

    auto gemm_mult0 = std::chrono::steady_clock::now();
    madness::cblas::gemm(madness::cblas::CBLAS_TRANSPOSE::NoTrans,
                         madness::cblas::CBLAS_TRANSPOSE::NoTrans, C.rows(),
                         C.cols(), A.cols(), 1.0, A.data(), A.rows(), B.data(),
                         B.rows(), 1.0, C.data(), C.rows());
    auto gemm_mult1 = std::chrono::steady_clock::now();
    double gemm_time
        = std::chrono::duration_cast
          <std::chrono::duration<double>>(gemm_mult1 - gemm_mult0).count();

    std::cout << input << ", " << percent_zero << ", " << dense_time << ", "
              << sparse_time << ", " << gemm_time << "\n";

    return 0;
}
