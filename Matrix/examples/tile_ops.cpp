#include "../tile_ops.h"
#include "../../tests/create_low_rank_array.h"
#include <ostream>

int main(int argc, char *argv[]) {
    int start = 50;
    int end = 501;

    std::string file_name = "FullRank_Square_gemm_no_add.dat";
    std::ofstream out_file(file_name);
    for (auto i = start; i < end; i += 50) {
        for (auto j = start; j < end; j += 50) {
            std::cout << "Full Gemm no add " << i << " " << j << std::endl;
            for (auto k = start; k < end; k += 50) {
                out_file << i << "," << j << "," << k << ",";
                FullRankTile<double> A{Eigen::MatrixXd::Random(i, k)};
                FullRankTile<double> B{Eigen::MatrixXd::Random(k, j)};
                double time = madness::wall_time();
                auto C = tile_ops::gemm(A, B, 1.0);
                time = madness::wall_time() - time;
                out_file << time << std::endl;
            }
        }
    }
    out_file.close();

    file_name = "FullRank_Square_gemm_with_add.dat";
    out_file.open(file_name);
    for (auto i = start; i < end; i += 50) {
        for (auto j = start; j < end; j += 50) {
            std::cout << "Full Gemm with add " << i << " " << j << std::endl;
            FullRankTile<double> C{Eigen::MatrixXd::Random(i, j)};
            for (auto k = start; k < end; k += 50) {
                out_file << i << "," << j << "," << k << ",";
                FullRankTile<double> A{Eigen::MatrixXd::Random(i, k)};
                FullRankTile<double> B{Eigen::MatrixXd::Random(k, j)};
                double time = madness::wall_time();
                tile_ops::gemm(C, A, B, 1.0, 1.0);
                time = madness::wall_time() - time;
                out_file << time << std::endl;
            }
        }
    }
    out_file.close();

    using matrix = typename LowRankTile<double>::template Matrix<double>;

    file_name = "LowRank_Square_gemm_no_add.dat";
    out_file.open(file_name);
    for (auto i = start; i < end; i += 50) {
        std::cout << "LR gemm no add " << i << " " << i << std::endl;
        for (auto r = 2; r < i / 2; r += 10) {
            out_file << i << "," << i << "," << i << "," << r << ",";
            LowRankTile<double> A{matrix::Random(i, r), matrix::Random(r, i)};
            LowRankTile<double> B{matrix::Random(i, r), matrix::Random(r, i)};
            double time = madness::wall_time();
            auto C = tile_ops::gemm(A, B, 1.0);
            time = madness::wall_time() - time;
            out_file << time << std::endl;
        }
    }
    out_file.close();

    file_name = "LowRank_Square_gemm_with_add.dat";
    out_file.open(file_name);
    for (auto i = start; i < end; i += 50) {
        std::cout << "LR gemm with add " << i << " " << i << std::endl;
        for (auto r = 2; r < i / 2; r += 10) {
            out_file << i << "," << i << "," << i << "," << r << ",";
            LowRankTile<double> A{matrix::Random(i, r), matrix::Random(r, i)};
            LowRankTile<double> B{matrix::Random(i, r), matrix::Random(r, i)};
            LowRankTile<double> C{matrix::Random(i, r), matrix::Random(r, i)};
            double time = madness::wall_time();
            tile_ops::gemm(C, A, B, 1.0, 1.0);
            time = madness::wall_time() - time;
            out_file << time << std::endl;
        }
    }
    out_file.close();

    file_name = "LowRank_Square_gemm_with_rounded_add.dat";
    out_file.open(file_name);
    for (auto i = start; i < end; i += 50) {
        for (auto j = start; j < end; j += 50) {
            std::cout << "LR gemm with rouned add " << i << " " << j
                      << std::endl;
            for (auto k = start; k <= std::min(i, j); k += 50) {
                auto max_rank = std::min({i, j, k}) / 2;
                for (auto r = 2; r < max_rank; r += 10) {
                    out_file << i << "," << j << "," << k << "," << r << ",";
                    LowRankTile
                        <double> A{matrix::Random(i, r), matrix::Random(r, k)};
                    LowRankTile
                        <double> B{matrix::Random(k, r), matrix::Random(r, j)};
                    LowRankTile
                        <double> C{matrix::Random(i, r), matrix::Random(r, j)};
                    double time = madness::wall_time();
                    tile_ops::gemm(C, A, B, 1.0, 1.0);
                    tile_ops::compress(C, 1e-07);
                    time = madness::wall_time() - time;
                    out_file << time << std::endl;
                }
            }
        }
    }
    out_file.close();

    return 0;
}
