
#include "../tile_algebra.h"

Eigen::MatrixXd scaled_rank_mat(int rows, int cols, double scale) {
    Eigen::MatrixXd M = Eigen::MatrixXd::Random(rows, cols);
    M *= 1000;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeThinU
                                             | Eigen::ComputeThinV);
    auto vals = svd.singularValues();
    auto iter = 1.0;
    for (auto i = vals.data(); i != vals.data() + vals.size(); ++i) {
        *i = *i / iter;
        iter *= scale;
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

double f_norm(Eigen::MatrixXd const &mat){
  auto squared_sum = 0.0;
  for(auto i = 0; i < mat.size(); ++i){
    auto elem = *(mat.data() + i);
    squared_sum += elem * elem;
  }
  return std::sqrt(squared_sum);
}


int main() {
    auto mat_size = 100;
    auto A = scaled_rank_mat(mat_size, mat_size, 2);
    auto B = scaled_rank_mat(mat_size, mat_size, 3);
    auto C = scaled_rank_mat(mat_size, mat_size, 4);
    Eigen::MatrixXd LA, RA, LB, RB, LC, RC;
    bool all_low = false;
    if (false && !algebra::Decompose_Matrix(A, LA, RA, 1e-07)
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

        std::cout << "F norm diff in A = " << diff_norm_A << std::endl;
        std::cout << "F norm diff in B = " << diff_norm_B << std::endl;
        std::cout << "F norm diff in C = " << diff_norm_C << std::endl;
        std::cout << "F norm diff in AB = " << (A * B - ABapprox).norm()
                  << " upper bound = " << AB_norm_upper << std::endl;
        C = A * B + C;
        QRCapprox = ABapprox + QRCapprox;
        std::cout << "F norm diff in (C = AB + C) = " << (C - QRCapprox).norm()
                  << " upper bound = " << C_norm_upper << std::endl;
    }

    Eigen::MatrixXd correct = scaled_rank_mat(100, 100, 2);
    Eigen::MatrixXd Lc, Rc;
    algebra::Decompose_Matrix(correct, Lc, Rc, 1e-07);
    Eigen::MatrixXd approx = Lc * Rc;
    auto norm_diff = f_norm(correct - approx);
    std::cout << "Initial diff = " << norm_diff << std::endl;
    auto iter = 1;
    auto cumulative_upper_bound = 0.0;
    while (norm_diff < 1e-03) {
        auto A = scaled_rank_mat(100, 100, 3);
        auto B = scaled_rank_mat(100, 100, 3);
        Eigen::MatrixXd LA, RA, LB, RB;
        algebra::Decompose_Matrix(A, LA, RA, 1e-07);
        algebra::Decompose_Matrix(B, LB, RB, 1e-07);

        Eigen::MatrixXd Aa = LA * RA;
        Eigen::MatrixXd Ba = LB * RB;
        auto A_diff = f_norm(A - Aa);
        auto B_diff = f_norm(B - Ba);

        std::cout << "\tA diff = " << A_diff << std::endl;
        std::cout << "\tB diff = " << B_diff << std::endl;

        auto upper_bound_mult = f_norm(A) * B_diff + f_norm(B) * A_diff;
        auto upper_bound_gemm = norm_diff + upper_bound_mult;
        cumulative_upper_bound += upper_bound_gemm;

        correct += A * B;
        approx += LA * RA * LB * RB;
        norm_diff = f_norm(correct - approx);
        std::cout << "On iteration " << iter << " diff = " << norm_diff
                  << " upper bound was " << upper_bound_gemm
                  << " cumulative upper bound was " << cumulative_upper_bound
                  << std::endl;
        ++iter;
    }


    return 0;
}
