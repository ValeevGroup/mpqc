#include "../include/tiledarray.h"

#include "../common/namespaces.h"
#include "../common/typedefs.h"

#include "../tensor/tcc_tile.h"

#include <iostream>

int main(int argc, char **argv){
    auto &world = madness::initialize(argc, argv);

    auto tr1 = TA::TiledRange1{2,2,2};
    auto tr = TA::TiledRange{tr1, tr1};

    TA::Array<double, 2, Tile<TA::Tensor<double>>> A(world, tr);
    A.set_all_local(1.0);

    std::cout << "Array A\n" << A << std::endl;

    madness::finalize();
    return 0;
}
