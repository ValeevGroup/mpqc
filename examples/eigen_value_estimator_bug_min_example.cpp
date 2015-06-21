#include <iostream>
#include <vector>

#include "../ta_routines/sqrt_inv.h"
#include "../ta_routines/array_to_eigen.h"
#include "../include/tiledarray.h"

#include "../utility/parallel_break_point.h"

namespace TA = TiledArray;

int main(int argc, char **argv) {
    auto &world = madness::initialize(argc, argv);

    if (argc < 4) {
        if (world.rank() == 0) {
            std::cout << "Program input is: ./<prog> <size> <#blocks> <repeats>"
                      << std::endl;
        }
        return 0;
    }

    auto size = std::stoi(argv[1]);
    auto nblocks = std::stoi(argv[2]);
    auto repeats = std::stoi(argv[3]);
    volatile int debug = (argc >= 5) ? std::stoi(argv[3]) : 0;

    tcc::utility::parallal_break_point(world,debug);

    auto block_size = size / nblocks;

    std::vector<int> blocking;
    for (auto i = 0; i < size; i += block_size) {
        blocking.push_back(i);
    }
    blocking.push_back(size);

    auto trange1 = TA::TiledRange1(blocking.begin(), blocking.end());
    auto trange = TA::TiledRange({trange1, trange1});

    auto shape_tensor = TA::Tensor<float>(trange.tiles(), 0.0);
    auto tsize = shape_tensor.range().size();
    auto tmap = TA::eigen_map(shape_tensor, tsize[0], tsize[1]);
    tmap(0, 0) = 1.0;
    for (auto i = 1u; i < tsize[0] - 1; ++i) {
        tmap(i, i) = 1.0;
        tmap(i - 1, i) = 0.5;
        tmap(i, i - 1) = 0.5;
    }
    tmap(tsize[0] - 1, tsize[0] - 1) = 1.0;


    if (world.rank() == 0) {
        std::cout << "Setting local tiles" << std::endl;
    }
    TA::Array<double, 2, TA::Tensor<double>, TA::SparsePolicy> array(
          world, trange, TA::SparseShape<float>(world, shape_tensor, trange));
    array.set_all_local(1.0);

    auto e_array = tcc::array_ops::array_to_eigen(array);
    if (world.rank() == 0) {
        std::cout << "Calculating correct eval range\n";
        Eigen::SelfAdjointEigenSolver<decltype(e_array)> es(e_array);
        std::cout << "Correct smallest eval is " << es.eigenvalues()[0]
                  << std::endl << " and smallest is "
                  << es.eigenvalues()[es.eigenvalues().size() - 1] << std::endl;
    }

    if (world.rank() == 0) {
        std::cout << "Starting Iterations" << std::endl;
    }
    for (auto i = 0; i < repeats; ++i) {
        auto evals = tcc::pure::eval_guess(array);
        if (world.rank() == 0) {
            std::cout << "Estimated evals on run " << i << " are " << evals[0]
                      << " and " << evals[1] << std::endl;
        }
    }

    world.gop.fence();
    madness::finalize();
    return 0;
}
