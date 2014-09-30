#include "../low_rank_tile.h"
#include "../tile_ops.h"
#include "../tile_mutations.h"
#include <iostream>

using matrix = Eigen::MatrixXd;

int main(int argc, char *argv[]) {
    auto dim = 0;
    auto rank = 0;
    if (argc >= 3) {
        dim = std::stoi(argv[1]);
        rank = std::stoi(argv[2]);
    } else {
        std::cout << "input should be ./program dim rank\n";
        return -1;
    }

    matrix A(dim,rank);
    matrix B(rank,dim);

    auto size = dim * rank;
    auto Ap = A.data();
    for(auto i = 0; i < size; ++i){
      *(Ap + i) = 1.0;
    }
    auto Bp = B.data();
    for(auto i = 0; i < size; ++i){
      *(Bp + i) = 1.0;
    }


    LowRankTile<double> lrA{A, B};
    LowRankTile<double> lrB{std::move(A), std::move(B)};

    double time = madness::wall_time();
    auto C = binary_ops::gemm(lrA, lrB, 1.0);
    time = madness::wall_time() - time;

    double time2 = madness::wall_time();
    ternary_mutations::gemm(C,lrA,lrB,1.0,1.0);
    time2 = madness::wall_time() - time2;

    std::cout << "Time for " << dim << "x" << dim
              << " gemm no add with rank = " << rank << " was " << time << " s\n"
              << "Time for gemm with add = " << time2 << " s\n";

    return 0;
}
