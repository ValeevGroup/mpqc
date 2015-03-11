#include "common/namespaces.h"
#include "common/typedefs.h"
#include "include/tiledarray.h"
#include "tensor/tile_pimpl.h"

int main(int argc, char *argv[]) {
    auto &world = madness::initialize(argc, argv);

    TA::TiledRange1 tr1(0,2,4);
    TA::TiledRange tr{tr1, tr1, tr1, tr1};
    TAArray<4, TA::DensePolicy> T{world, std::move(tr)};
    T.fill_local(1.0);

    Eig::MatrixXd M = Eig::MatrixXd::Random(4,4);
    
    std::cout << "Tensor = \n" << T << "\n" << std::endl;
    std::cout << "Matrix = \n" << M << std::endl;
    return 0;
}
