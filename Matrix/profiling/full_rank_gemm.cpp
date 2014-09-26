#include "../full_rank_tile.h"
#include "../tile_ops.h"
#include "../tile_mutations.h"
#include "../include/tiledarray.h"
#include <iostream>
#include <string>

using matrix = Eigen::MatrixXd;

int main(int argc, char *argv[]) {
    auto dim = 0;
    if (argc >= 2) {
        dim = std::stoi(argv[1]);
    } else {
        std::cout << "input should be ./program dim\n";
        return -1;
    }

    matrix A(dim, dim);
    TiledArray::Tensor<double> TA_A(TiledArray::Range(dim, dim), 1.0);

    TiledArray::math::GemmHelper gemm_helper(
        madness::cblas::CBLAS_TRANSPOSE::NoTrans,
        madness::cblas::CBLAS_TRANSPOSE::NoTrans, 2, 2, 2);


    auto size = dim * dim;
    auto Ap = A.data();
    for (auto i = 0; i < size; ++i) {
        *(Ap + i) = 1.0;
    }

    FullRankTile<double> fA{A.eval()};
    FullRankTile<double> fB{std::move(A)};

    double time = madness::wall_time();
    auto C = binary_ops::gemm(fA, fB, 1.0);
    time = madness::wall_time() - time;

    double ta_time = madness::wall_time();
    TA_A = TA_A.gemm(TA_A, 1.0, gemm_helper);
    ta_time = madness::wall_time() - ta_time;

    double time2 = madness::wall_time();
    ternary_mutations::gemm(C, fA, fB, 1.0, 1.0);
    time2 = madness::wall_time() - time2;

    std::cout << "Time for " << dim << "x" << dim << " gemm no add was " << time
              << " s\n"
              << "Time for gemm with add = " << time2 << " s\n"
              << "Time for TA gemm no add = " << ta_time << " s\n";

    return 0;
}
