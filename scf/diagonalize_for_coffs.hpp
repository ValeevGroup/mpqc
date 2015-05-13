#pragma once
#ifndef TCC_SCF_DIAGONALIZEFORCOFFS_H
#define TCC_SCF_DIAGONALIZEFORCOFFS_H

#include "../common/namespaces.h"

#include "../include/tiledarray.h"
#include "../include/eigen.h"

#include "../array_ops/array_to_eigen.h"

namespace tcc {
namespace scf {

using Array2 = TA::Array<double, 2, TA::Tensor<double>, TA::SparsePolicy>;

TA::TiledRange1 tr_occupied(int guess, int occ) {
    auto nblocks = (guess < occ) ? guess : occ;
    auto block_size = std::max(occ / nblocks, 1);
    std::vector<std::size_t> blocks;
    blocks.reserve(nblocks + 1);
    blocks.push_back(0);
    for (auto i = block_size; i < occ; i += block_size) {
        blocks.push_back(i);
    }
    blocks.push_back(occ);
    return TA::TiledRange1(blocks.begin(), blocks.end());
}

Array2 Coeffs_from_fock(Array2 const &F, Array2 const &S, TA::TiledRange1 tr_i,
                        unsigned int occ) {
    auto F_eig = array_ops::array_to_eigen(F);
    auto S_eig = array_ops::array_to_eigen(S);

    Eig::GeneralizedSelfAdjointEigenSolver<decltype(S_eig)> es(F_eig, S_eig);
    decltype(S_eig) C = es.eigenvectors().leftCols(occ);

    auto tr_ao = S.trange().data()[0];

    return array_ops::eigen_to_array<TA::Tensor<double>>(S.get_world(), C,
                                                         tr_ao, tr_i);
}

} // namespace scf
} // namespace tcc


#endif // TCC_SCF_DIAGONALIZEFORCOFFS_H
