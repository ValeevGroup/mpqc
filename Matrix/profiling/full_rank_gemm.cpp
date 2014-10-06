#include "../full_rank_tile.h"
#include "../tile_ops.h"
#include "../tile_mutations.h"
#include <iostream>
#include <string>

using matrix = Eigen::MatrixXd;

int main(int argc, char *argv[]) {
    auto dim = 0;
    auto rank = 0;
    if (argc >= 2) {
        dim = std::stoi(argv[1]);
    } else {
        std::cout << "input should be ./program dim\n";
        return -1;
    }

    matrix A(dim,dim);

    auto size = dim * rank;
    auto Ap = A.data();
    for(auto i = 0; i < size; ++i){
      *(Ap + i) = 1.0;
    }

    FullRankTile<double> fA{A};
    FullRankTile<double> fB{std::move(A)};

    double time = madness::wall_time();
    auto C = binary_ops::gemm(fA, fB, 1.0);
    time = madness::wall_time() - time;

    double time2 = madness::wall_time();
    C = ternary_mutations::gemm(std::move(C),fA,fB,1.0,1.0);
    time2 = madness::wall_time() - time2;

    std::cout << "Time for " << dim << "x" << dim
              << " gemm no add was " << time << " s\n"
              << "Time for gemm with add = " << time2 << " s\n";

    return 0;
}
