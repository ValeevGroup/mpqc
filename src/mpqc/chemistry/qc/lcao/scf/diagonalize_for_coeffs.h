
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_DIAGONALIZE_FOR_COEFFS_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_DIAGONALIZE_FOR_COEFFS_H_


#include "mpqc/math/tensor/clr/array_to_eigen.h"

#include <tiledarray.h>
#include "mpqc/math/external/eigen/eigen.h"

#include "mpqc/chemistry/qc/lcao/integrals/make_engine.h"

#include "mpqc/math/tensor/clr/vector_localization.h"
#include "mpqc/math/tensor/clr/decomposed_tensor_algebra.h"


namespace mpqc {
namespace lcao {
namespace scf {

using Array2 = TA::Array<double, 2, TA::Tensor<double>, TA::SparsePolicy>;

inline TA::TiledRange1 tr_occupied(int guess, int occ) {
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

inline Array2 Coeffs_from_fock(Array2 const &F, Array2 const &S, TA::TiledRange1 tr_i,
                        unsigned int occ, unsigned int nocc_clusters, 
                        bool use_chol_vectors = false) {
    auto F_eig = math::array_to_eigen(F);
    auto S_eig = math::array_to_eigen(S);

    Eigen::GeneralizedSelfAdjointEigenSolver<decltype(S_eig)> es(F_eig, S_eig);
    decltype(S_eig) C = es.eigenvectors().leftCols(occ);

    if (use_chol_vectors) {
        decltype(S_eig) D = C * C.transpose();
        unsigned int rank = math::piv_cholesky(D);
        if (rank == occ) {
            C = D;
        } else {
            std::cout << "Cholesky Rank was " << rank << std::endl;
            throw std::runtime_error("The rank of the cholesky vectors was not "
                                     "equal to the number of occupied "
                                     "orbitals");
        }
    }

    auto tr_ao = S.trange().data()[0];

    return math::eigen_to_array<TA::Tensor<double>,TA::SparsePolicy>(S.world(), C,
                                                         tr_ao, tr_i);
}

}  // namespace  scf
}  // namespace  lcao
}  // namespace  mpqc


#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_DIAGONALIZE_FOR_COEFFS_H_
