#include "../include/tiledarray.h"

#include "../common/namespaces.h"
#include "../common/typedefs.h"

#include "../utility/parallel_break_point.h"

#include "../tensor/tcc_tile.h"

#include <iostream>

int main(int argc, char **argv) {
    auto &world = madness::initialize(argc, argv);
    auto debug = (argc >= 2) ? std::stoi(argv[1]) : 0;

    tcc::utility::parallal_break_point(world, debug);

    auto tr1 = TA::TiledRange1{0, 2, 4};
    auto tr = TA::TiledRange{tr1, tr1};

    TA::Array<double, 2, Tile<TA::Tensor<double>>, TA::DensePolicy> A(world,
                                                                      tr);
    A.set_all_local(3.0);

    std::cout << "Array A\n" << A << std::endl;

    TA::Array<double, 2, Tile<TA::Tensor<double>>, TA::DensePolicy> B(world,
                                                                      tr);
    B.set_all_local(2.0);

    std::cout << "B = \n" << B << std::endl;


    decltype(A) C;
    C("i,j") = A("i,j") + B("i,j");

    std::cout << "C = \n" << C << std::endl;

    C("i,j") = C("i,j") + A("i,j");
    std::cout << "C = \n" << C << std::endl;

    C("i,j") = C("i,j") + C("j,i");
    std::cout << "C = \n" << C << std::endl;

    C("j,i") = B("i,j") + A("i,j");
    std::cout << "C = \n" << C << std::endl;

    C("i,j") = B("i,k") * A("k,j");
    std::cout << "C = \n" << C << std::endl;

    madness::finalize();
    return 0;
}
