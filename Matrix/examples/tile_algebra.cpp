#include "../tile_algebra.h"
#include "../../tests/create_low_rank_array.h"
#include <ostream>

int main(int argc, char *argv[]) {
    int start = 50;
    int end = 2001;

    Eigen::MatrixXd A;
    Eigen::MatrixXd B;
    Eigen::MatrixXd C;

    for (auto i = start; i < end; i += 50) {
        for (auto j = start; j < end; j += 50) {
            std::string file_name = "Square_gemm_no_add_" + std::to_string(i)
                                    + "_" + std::to_string(j) + ".dat";
            std::ofstream out_file(file_name);
            for (auto k = start; k < end; k += 50) {
                out_file << i << "," << j << "," << k << ",";
                double time = madness::wall_time();
                A = Eigen::MatrixXd::Random(i, k);
                B = Eigen::MatrixXd::Random(k, j);
                C = algebra::cblas_gemm(A, B, 1.0);
                time = madness::wall_time() - time;
                out_file << time << std::endl;
            }
            out_file.close();
            std::cout << "Done with Gemm no add " << i << " " << j << std::endl;
        }
    }

    // C square Gemm
    for (auto i = start; i < end; i += 50) {
        for (auto j = start; j < end; j += 50) {
            std::string file_name = "Square_gemm_with_add_" + std::to_string(i)
                                    + "_" + std::to_string(j) + ".dat";
            std::ofstream out_file(file_name);
            C = Eigen::MatrixXd::Random(i, j);
            for (auto k = start; k < end; k += 50) {
                out_file << i << "," << j << "," << k << ",";
                double time = madness::wall_time();
                A = Eigen::MatrixXd::Random(i, k);
                B = Eigen::MatrixXd::Random(k, j);
                algebra::cblas_gemm_inplace(A, B, C, 1.0, 1.0);
                time = madness::wall_time() - time;
                out_file << time << std::endl;
            }
            out_file.close();
            std::cout << "Done with Gemm " << i << " " << j << std::endl;
        }
    }

    for (auto i = start; i < end; i += 50) {
        for (auto j = start; j < end; j += 50) {
            std::string file_name = "QR_Init_" + std::to_string(i) + "_"
                                    + std::to_string(j) + ".dat";
            std::ofstream out_file(file_name);
            for (auto k = 10; k < std::min(i, j); k += 10) {
                C = HeirChemTest::low_rank_matrix<double>(i, j, k);
                out_file << i << "," << j << "," << k << ",";
                double time = madness::wall_time();
                algebra::QrInit(C, A, B, 1e-09);
                time = madness::wall_time() - time;
                out_file << time << std::endl;
            }
            out_file.close();
            std::cout << "Done with Qr " << i << " " << j << std::endl;
        }
    }

    return 0;
}
