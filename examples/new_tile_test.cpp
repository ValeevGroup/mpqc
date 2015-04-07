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

    // Make low rank tensor
    TA::Range r{50, 10, 10};
    TA::Tensor<double> M(r);
    auto L = TA::EigenMatrixXd::Random(50, 10);
    auto R = TA::EigenMatrixXd::Random(10, 100);
    auto map_M = TA::eigen_map(M, 50, 100);
    map_M = L * R;
    tensor::DecomposedTensor<double> lr_tensor{1e-7, std::move(M)};

    // Decompose it with QR
    auto output_T = tensor::algebra::two_way_decomposition(lr_tensor);

    if (output_T.empty()) {
        std::cout << "Something be broken" << std::endl;
    } else {
        std::cout << "Maybe this actually worked" << std::endl;
        std::cout << "Output rank = " << output_T.rank() << std::endl;

        auto out = tensor::algebra::combine(output_T);

        auto diff = lr_tensor.tensor(0).subt(out).norm();
        std::cout << "The norm diff of the tensors was " << diff << std::endl;
    }

    auto tr1 = TA::TiledRange1{0, 2, 4};
    auto tr = TA::TiledRange{tr1, tr1};

    TA::Array<double, 2, Tile<tensor::DecomposedTensor<double>>,
              TA::DensePolicy> A(world, tr);
    A.set_all_local(1.0);

    std::cout << A << std::endl;


    /* decltype(A) B; */
    /* B("i,j") = A("i,j") + A("i,j"); */

    /* decltype(A) C; */
    /* C("i,j") = A("i,k") + B("k,j"); */

    madness::finalize();
    return 0;
}
