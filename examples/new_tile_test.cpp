#include "../include/tiledarray.h"

#include "../common/namespaces.h"
#include "../common/typedefs.h"

#include "../tensor/tcc_tile.h"

#include <iostream>

int main(int argc, char **argv) {
    auto &world = madness::initialize(argc, argv);

    auto tr1 = TA::TiledRange1{0, 2, 4};
    auto tr = TA::TiledRange{tr1, tr1};

    TA::Array<double, 2, Tile<TA::Tensor<double>>, TA::DensePolicy> A(world,
                                                                      tr);

    A.set_all_local(1.0);

    std::cout << "Array A\n" << A << std::endl;

    for (auto fut : A) {
        std::cout << "Tile is empty? " << fut.get().empty() << "\n";
        std::cout << "Norm of tile = " << fut.get().norm() << "\n";
    }

    TA::Array<double, 2, Tile<TA::Tensor<double>>, TA::DensePolicy> B(world,
                                                                      tr);
    B.set_all_local(2.0);

    decltype(A) C;
    C("i,j") = A("i,j") + B("i,j");

    std::cout << "C = \n" << C << std::endl;

    C("i,j") = C("i,j") + A("i,j");

    std::cout << "C = \n" << C << std::endl;

    C("i,j") = C("i,j") + C("j,i");

    std::cout << "C = \n" << C << std::endl;

    madness::finalize();
    return 0;
}
