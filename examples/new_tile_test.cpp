#include "../include/tiledarray.h"

#include "../common/namespaces.h"
#include "../common/typedefs.h"

#include "../utility/parallel_break_point.h"

#include "../tensor/decomposed_tensor_nonintrusive_interface.h"
#include "../tensor/tcc_tile.h"

#include <iostream>

using namespace tcc;
int main(int argc, char **argv) {
    auto &world = madness::initialize(argc, argv);
    auto debug = (argc >= 2) ? std::stoi(argv[1]) : 0;

    tcc::utility::parallal_break_point(world, debug);

    auto tr1 = TA::TiledRange1{0, 2, 4};
    auto tr = TA::TiledRange{tr1, tr1};

    TA::Array<double, 2, Tile<tensor::DecomposedTensor<double>>,
              TA::DensePolicy> A(world, tr);


    decltype(A) B;
    B("i,j") = A("i,j") + A("i,j");

    decltype(A) C;
    C("i,j") = A("i,k") + B("k,j");

    madness::finalize();
    return 0;
}
