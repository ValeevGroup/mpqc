#include "../include/tiledarray.h"

#include "../common/namespaces.h"
#include "../common/typedefs.h"

#include "../utility/parallel_break_point.h"

#include "../tensor/decomposed_tensor_nonintrusive_interface.h"
#include "../tensor/tcc_tile.h"
#include "../tensor/density_tensor.h"
#include "../utility/time.h"

#include <iostream>

using namespace tcc;

TA::Tensor<double> lr_ta_tensor(std::size_t dim1, std::size_t dim2,
                                std::size_t dim3, std::size_t rank) {
    TA::Range r{dim1, dim2, dim3};
    TA::Tensor<double> M(r);
    auto L = TA::EigenMatrixXd::Random(dim1, rank);
    auto R = TA::EigenMatrixXd::Random(rank, dim2 * dim3);
    auto map_M = TA::eigen_map(M, dim1, dim2 * dim3);
    map_M = L * R;

    return M;
}

TA::Tensor<double> random_matrix(std::size_t dim) {
    TA::Range r{dim, dim};
    TA::Tensor<double> M(r);
    auto map_M = TA::eigen_map(M, dim, dim);
    map_M = TA::EigenMatrixXd::Random(dim, dim);

    return M;
}

void
do_test(tcc::tensor::DecomposedTensor<double> &lr,
        tcc::tensor::DecomposedTensor<double> const &d,
        tcc::tensor::DecomposedTensor<double> &output,
        TA::math::GemmHelper const &gh, TA::Tensor<double> const &correct) {
    if (output.empty()) {
        auto lr_copy = tensor::clone(lr);
        auto total_time = 0.0;
        auto t0 = tcc_time::now();
        tensor::gemm(lr, lr_copy, d, 1.0, gh);
        auto t1 = tcc_time::now();
        total_time += tcc_time::duration_in_s(t0, t1);
        auto my_product = tensor::algebra::combine(lr);
        auto diff = my_product.subt(correct).norm();
        std::cout << "No Decomp diff = " << diff << std::endl;
        std::cout << "Time me = " << total_time << std::endl;
    } else {
        auto o_copy = tensor::clone(output);
        auto total_time = 0.0;
        auto t0 = tcc_time::now();
        tensor::gemm(output, o_copy, d, 1.0, gh);
        auto t1 = tcc_time::now();
        total_time += tcc_time::duration_in_s(t0, t1);
        auto my_product = tensor::algebra::combine(output);
        auto diff = my_product.subt(correct).norm();
        std::cout << "Decomp diff = " << diff << std::endl;
        std::cout << "Time me = " << total_time << std::endl;
    }
}

int main(int argc, char **argv) {
    auto &world = madness::initialize(argc, argv);
    auto df_dim = std::stoi(argv[1]);
    auto bs_dim = std::stoi(argv[2]);
    auto rank = std::stoi(argv[3]);

    auto TA_lr = lr_ta_tensor(df_dim, bs_dim, bs_dim, rank);
    auto D_TA = random_matrix(bs_dim);
    tensor::DecomposedTensor<double> lr_tensor{1e-7, TA_lr};
    tensor::DecomposedTensor<double> d_tensor{1e-7, D_TA};
    auto output_t = tensor::algebra::two_way_decomposition(lr_tensor);

    const auto NoT = madness::cblas::CBLAS_TRANSPOSE::NoTrans;
    auto gh = TA::math::GemmHelper(NoT, NoT, 3, 3, 2);

    // dummy call
    TA_lr.gemm(D_TA, 1.0, gh);

    auto ta_lr_copy = TA_lr.clone();
    auto t0 = tcc_time::now();
    TA_lr.gemm(ta_lr_copy, D_TA, 1.0, gh);
    auto t1 = tcc_time::now();
    auto total_time = tcc_time::duration_in_s(t0, t1);
    std::cout << "Ta time = " << total_time << std::endl;

    do_test(lr_tensor, d_tensor, output_t, gh, TA_lr);

    madness::finalize();
    return 0;
}
