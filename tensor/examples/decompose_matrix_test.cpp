
#include "../tile_algebra.h"
#include <ctime>
#include <iomanip>

Eigen::MatrixXd scaled_rank_mat(int rows, int cols, double scale) {
    Eigen::MatrixXd M = Eigen::MatrixXd::Random(rows, cols);
    for (auto i = 0; i < M.rows(); ++i) {
        for (auto j = 0; j < M.cols(); ++j) {
            M(i, j) = std::abs(M(i, j));
        }
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeThinU
                                             | Eigen::ComputeThinV);
    auto vals = svd.singularValues();
    auto iter = 1.0;
    auto counter = 1;
    for (auto i = vals.data(); i != vals.data() + vals.size(); ++i) {
        *i = *i / iter;
        iter = std::pow(counter, scale); // scale;
        ++counter;
    }
    return svd.matrixU() * vals.asDiagonal() * svd.matrixV().transpose();
}


Eigen::MatrixXd svd_approx_mat(Eigen::MatrixXd const &mat, double cut) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(mat, Eigen::ComputeThinU
                                               | Eigen::ComputeThinV);
    auto vals = svd.singularValues();
    for (auto i = vals.data(); i != vals.data() + vals.size(); ++i) {
        if (*i < cut) {
            *i = 0;
        }
    }
    return svd.matrixU() * vals.asDiagonal() * svd.matrixV().transpose();
}

double f_norm(Eigen::MatrixXd const &mat) {
    auto squared_sum = 0.0;
    for (auto i = 0; i < mat.size(); ++i) {
        auto elem = *(mat.data() + i);
        squared_sum += elem * elem;
    }
    return std::sqrt(squared_sum);
}


int main(int argc, char **argv) {
    srand(time(NULL));
    auto mat_size = 0;
    auto scaling_factor = 5.0;
    if (argc < 2) {
        std::cout << "Usage is ./prog matrix_size scaling_factor(optional)\n";
        return 1;
    }
    if (argc >= 2) {
        mat_size = std::stoi(argv[1]);
        if (mat_size <= 0) {
            std::cout << "matrix_size must be a posative int\n";
            return 1;
        }
    }
    if (argc >= 3) {
        scaling_factor = std::stod(argv[2]);
        if (scaling_factor <= 0) {
            std::cout << "scaling_factor must be a posative floating type\n";
            return 1;
        }
    }
    std::cout << "Matrix size is " << mat_size << " scaling factor is "
              << scaling_factor << std::endl;
    auto A = scaled_rank_mat(mat_size, mat_size, scaling_factor);
    auto B = scaled_rank_mat(mat_size, mat_size, scaling_factor);
    auto C = scaled_rank_mat(mat_size, mat_size, scaling_factor);
    Eigen::MatrixXd LA, RA, LB, RB, LC, RC;
    bool all_low = false;
    if (!algebra::Decompose_Matrix(A, LA, RA, 1e-07)
        && !algebra::Decompose_Matrix(B, LB, RB, 1e-07)
        && !algebra::Decompose_Matrix(C, LC, RC, 1e-07)) {
        all_low = true;
    }
    if (all_low) {
        Eigen::MatrixXd QRAapprox = LA * RA;
        Eigen::MatrixXd QRBapprox = LB * RB;
        Eigen::MatrixXd QRCapprox = LC * RC;
        Eigen::MatrixXd ABapprox = QRAapprox * QRBapprox;

        auto diff_norm_A = (A - QRAapprox).norm();
        auto diff_norm_B = (B - QRBapprox).norm();
        auto diff_norm_C = (C - QRCapprox).norm();


        auto norm_A = A.norm();
        auto norm_B = B.norm();

        auto AB_norm_upper = (norm_A * diff_norm_B) + norm_B * diff_norm_A;
        auto C_norm_upper = diff_norm_C + AB_norm_upper;

        std::cout << "Matrix dimensions are " << mat_size << "x" << mat_size
                  << " with threshold of " << 1e-07 << std::endl;
        std::cout << "F norm error in A = " << diff_norm_A << std::endl;
        std::cout << "F norm error B = " << diff_norm_B << std::endl;
        std::cout << "F norm error C = " << diff_norm_C << std::endl;
        std::cout << "F norm error AB = " << (A * B - ABapprox).norm()
                  << " upper bound for AB = " << AB_norm_upper << std::endl;
        C = A * B + C;
        QRCapprox = ABapprox + QRCapprox;
        std::cout << "F norm error in (C = AB + C) = " << (C - QRCapprox).norm()
                  << " upper bound = " << C_norm_upper << std::endl;
    }

    std::cout << "\n\nBegining iterative gemm to see how error accumulates "
                 "over many iterations of C = AB + C.\nA new A and B will be "
                 "created for each iteration.\nAll error is using F norms.\n"
                 "Stoping condition is when error >= 1e-05.\n\n";
    Eigen::MatrixXd correct
        = scaled_rank_mat(mat_size, mat_size, scaling_factor);
    Eigen::MatrixXd Lc, Rc;
    algebra::Decompose_Matrix(correct, Lc, Rc, 1e-07);
    Eigen::MatrixXd approx = Lc * Rc;
    auto norm_diff = f_norm(correct - approx);
    std::cout
        << "Initial error in C using a threshold of 1e-07 for matrices of dim "
        << mat_size << "x" << mat_size << " was " << norm_diff << std::endl;
    auto iter = 1;
    auto cumulative_upper_bound = 0.0;
    std::vector<double> iteration_error;
    while (norm_diff < 1e-05) {
        auto A = scaled_rank_mat(mat_size, mat_size, scaling_factor);
        auto B = scaled_rank_mat(mat_size, mat_size, scaling_factor);
        Eigen::MatrixXd LA, RA, LB, RB;
        algebra::Decompose_Matrix(A, LA, RA, 1e-07);
        algebra::Decompose_Matrix(B, LB, RB, 1e-07);

        Eigen::MatrixXd Aa = LA * RA;
        Eigen::MatrixXd Ba = LB * RB;
        auto A_diff = f_norm(A - Aa);
        auto B_diff = f_norm(B - Ba);

        std::cout << "\tA error = " << A_diff << std::endl;
        std::cout << "\tB error = " << B_diff << std::endl;

        auto upper_bound_mult = f_norm(A) * B_diff + f_norm(B) * A_diff;
        auto upper_bound_gemm = norm_diff + upper_bound_mult;
        cumulative_upper_bound += upper_bound_gemm;

        correct += A * B;
        approx += LA * RA * LB * RB;
        norm_diff = f_norm(correct - approx);
        iteration_error.emplace_back(norm_diff);
        std::cout << "On iteration " << iter << " acutal error = " << norm_diff
                  << " iteration upper bound was " << upper_bound_gemm
                  << " cumulative upper bound was " << cumulative_upper_bound
                  << std::endl;
        ++iter;
    }

    std::cout << "\nPrinting the actual error in each generation.\n";
    auto counter = 1;
    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout.setf(std::ios::showpoint);
    std::cout.precision(14);
    for (auto elem = iteration_error.begin(); elem != iteration_error.end();
         ++elem) {
        if (counter >= 5) {
            std::cout << "\n";
            counter = 1;
        }
        std::cout << *elem;
        if (elem != iteration_error.end() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "\n";


    return 0;
}
