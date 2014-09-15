#include "../tile_algebra.h"
#include "../../tests/create_low_rank_array.h"
#include <ostream>

int main(int argc, char *argv[]) {
    int start = 50;
    int end = 501;

    Eigen::MatrixXd A;
    Eigen::MatrixXd B;
    Eigen::MatrixXd C;

    std::string file_name = "Square_gemm_no_add.dat";
    std::ofstream out_file(file_name);
    for (auto i = start; i < end; i += 50) {
        for (auto j = start; j < end; j += 50) {
            std::cout << "Gemm no add " << i << " " << j << std::endl;
            for (auto k = start; k < end; k += 50) {
                out_file << i << "," << j << "," << k << ",";
                A = Eigen::MatrixXd::Random(i, k);
                B = Eigen::MatrixXd::Random(k, j);
                double time = madness::wall_time();
                C = algebra::cblas_gemm(A, B, 1.0);
                time = madness::wall_time() - time;
                out_file << time << std::endl;
            }
        }
    }
    out_file.close();

    file_name = "Square_gemm_with_add.dat";
    out_file.open(file_name);
    for (auto i = start; i < end; i += 50) {
        for (auto j = start; j < end; j += 50) {
            std::cout << "Gemm with add " << i << " " << j << std::endl;
            C = Eigen::MatrixXd::Random(i, j);
            for (auto k = start; k < end; k += 50) {
                out_file << i << "," << j << "," << k << ",";
                A = Eigen::MatrixXd::Random(i, k);
                B = Eigen::MatrixXd::Random(k, j);
                double time = madness::wall_time();
                algebra::cblas_gemm_inplace(A, B, C, 1.0, 1.0);
                time = madness::wall_time() - time;
                out_file << time << std::endl;
            }
        }
    }
    out_file.close();


    file_name = "ColPivotedQr.dat";
    out_file.open(file_name);
    for (auto i = start; i < end; i += 50) {
        for (auto j = start; j < end; j += 50) {
          std::cout << "Qr " << i << " " << j << std::endl;
            for (auto k = 10; k < std::min(i, j); k += 20) {
                C = TCC::test::low_rank_matrix<double>(i, j, k);
                out_file << i << "," << j << "," << k << ",";
                double time = madness::wall_time();
                algebra::ColPivotedQr(C, A, B, 1e-09);
                time = madness::wall_time() - time;
                out_file << time << std::endl;
            }
        }
    }
    out_file.close();

    return 0;
}
