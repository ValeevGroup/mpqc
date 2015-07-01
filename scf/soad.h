#pragma once
#ifndef TCC_INTEGRALS_SCF_SOAD_H
#define TCC_INTEGRALS_SCF_SOAD_H

#include "../basis/basis.h"
#include "../basis/basis_set.h"
#include "../molecule/cluster.h"
#include "../molecule/cluster_collapse.h"

#include "../include/tiledarray.h"
#include "../common/typedefs.h"
#include "../utility/make_array.h"
#include "../tensor/tcc_tile.h"
#include "../tensor/decomposed_tensor.h"
#include "../tensor/decomposed_tensor_nonintrusive_interface.h"

#include "../integrals/sparse_task_integrals.h"

#include "../ta_routines/array_to_eigen.h"

#include <memory>
#include <vector>

namespace tcc {
namespace integrals {
namespace scf {

using DataType = tensor::DecomposedTensor<double>;
using TileType = tensor::Tile<DataType>;
using ArrayType = TA::Array<double, 2, TileType, TA::SparsePolicy>;
using Eri3ArrayType = TA::Array<double, 3, TileType, TA::SparsePolicy>;


/*! \brief constructs a tile using Edwards SOAD strategy.
 *
 * This function creates a tile to initialize a density guess based on
 * the Superposition of Atomic Densities (SOAD) approach.  The function
 * expects to receive symmetric diagonal ranges.
 */
TileType soad_tile(std::shared_ptr<molecule::Cluster> cluster, TA::Range range,
                   double cut) {
    // Grab the range extent since it will be reused several times.
    const auto extent = range.extent();

    TA::Tensor<double> tensor(TA::Range{extent[0], extent[1]}, 0.0);
    auto t_map = TA::eigen_map(tensor, extent[0], extent[1]);

    auto atoms = molecule::collapse_to_atoms(*cluster);

    size_t ao_offset = 0; // first AO of this atom
    for (auto const &atom : atoms) {
        const auto Z = atom.charge();

        if (Z == 1 || Z == 2) {              // H, He
            t_map(ao_offset, ao_offset) = Z; // all electrons go to the 1s
            ao_offset += 1;
        } else if (Z <= 10) {
            t_map(ao_offset, ao_offset) = 2; // 2 electrons go to the 1s
            t_map(ao_offset + 1, ao_offset + 1) = (Z == 3) ? 1 : 2;
            // smear the remaining electrons in 2p orbitals
            const double num_electrons_per_2p = (Z > 4) ? (double)(Z - 4) / 3
                                                        : 0;
            for (auto xyz = 0; xyz != 3; ++xyz)
                t_map(ao_offset + 2 + xyz, ao_offset + 2 + xyz)
                      = num_electrons_per_2p;
            ao_offset += 5;
        }
    }

    // Multiply by 1/2 because we are using the convention that the Density is
    // idempotent.
    t_map *= 0.5;
    return TileType(range, DataType(cut, std::move(tensor)));
}


/*! \brief creates an initial guess for the density matrix using SOAD.
 *
 * This function creates a guess in the 3-21g basis using the Superposition
 * of atomic densities (SOAD) approach.
 */
ArrayType minimal_density_guess(
      madness::World &world,
      std::vector<std::shared_ptr<molecule::Cluster>> const &clusters,
      basis::Basis const &min_bs, double cut) {

    TA::TiledRange trange
          = sparse::create_trange(utility::make_array(min_bs, min_bs));

    // Make a shape tensor and set it such that only diagonal elements will 
    // be significant
    TA::Tensor<float> tile_norms(trange.tiles(), 0.0);
    auto const extent = tile_norms.range().extent();
    auto tn_map = TA::eigen_map(tile_norms, extent[0], extent[1]);
    for (auto i = 0; i < tn_map.rows(); ++i) {
        tn_map(i, i) = std::numeric_limits<float>::max(); // Go big 
    }

    // Create the array and shape
    TA::SparseShape<float> shape(world, tile_norms, trange);
    ArrayType D_min(world, trange, shape);

    // Initialize D_min tiles by calling soad_tile
    auto const &pmap = *D_min.get_pmap();
    auto it = pmap.begin();
    const auto end = pmap.end();
    for (; it != end; ++it) {
        if (!D_min.is_zero(*it)) {
            auto cluster_ord = trange.tiles().idx(*it)[0];
            madness::Future<tensor::Tile<tensor::DecomposedTensor<double>>> tile
                  = world.taskq.add(soad_tile, clusters[cluster_ord],
                                    trange.make_tile_range(*it), cut);
            D_min.set(*it, tile);
        }
    }

    D_min.truncate(); // Fix tile norms.
    return D_min;
}

template <typename SharedEnginePool, typename Op>
ArrayType fock_from_minimal(
      madness::World &world, basis::Basis const &obs, basis::Basis const &df_bs,
      SharedEnginePool eng_pool, ArrayType const &V_inv, ArrayType const &H,
      Eri3ArrayType const &W,
      std::vector<std::shared_ptr<molecule::Cluster>> const &clusters,
      double cut, Op op) {

    basis::BasisSet min_bs_set("sto-3g");

    std::streambuf *cout_sbuf = std::cout.rdbuf(); // Silence libint printing.
    std::ofstream fout("/dev/null");
    std::cout.rdbuf(fout.rdbuf());
    basis::Basis min_bs{min_bs_set.create_soad_basis(clusters)};
    std::cout.rdbuf(cout_sbuf);
    auto D_min = minimal_density_guess(world, clusters, min_bs, cut);

    auto EriJ = BlockSparseIntegrals(
          world, eng_pool, utility::make_array(df_bs, min_bs, min_bs), op);
    auto EriK
          = BlockSparseIntegrals(world, eng_pool,
                                 utility::make_array(df_bs, obs, min_bs), op);

    decltype(D_min) J, K, F;
    J("i,j") = W("X,i,j") * (EriJ("X,a,b") * D_min("a,b"));
    Eri3ArrayType W_K;
    W_K("X,a,i") = V_inv("X,P") * EriK("P,i,a");
    K("i,j") = W_K("X,a,i") * (EriK("X, j, b") * D_min("b,a"));
    F("i,j") = H("i,j") + 2 * J("i,j") - K("i,j");

    return F;
}

template <typename SharedEnginePool, typename Op>
ArrayType fock_from_minimal_v_oh(
      madness::World &world, basis::Basis const &obs, basis::Basis const &df_bs,
      SharedEnginePool eng_pool, ArrayType const &H, ArrayType const &V_inv_oh,
      Eri3ArrayType const &Xab,
      std::vector<std::shared_ptr<molecule::Cluster>> const &clusters,
      double cut, Op op) {

    basis::BasisSet min_bs_set("sto-3g");

    std::streambuf *cout_sbuf = std::cout.rdbuf(); // Silence libint printing.
    std::ofstream fout("/dev/null");
    std::cout.rdbuf(fout.rdbuf());
    basis::Basis min_bs{min_bs_set.create_soad_basis(clusters)};
    std::cout.rdbuf(cout_sbuf);
    auto D_min = minimal_density_guess(world, clusters, min_bs, cut);


    decltype(D_min) J, K, F, Coeffs;

    {
        auto EriJ = BlockSparseIntegrals(
              world, eng_pool, utility::make_array(df_bs, min_bs, min_bs), op);
        TA::Array<double, 1, typename decltype(EriJ)::value_type,
                  TA::SparsePolicy> trans;
        trans("P") = EriJ("P,a,b") * D_min("a,b");
        J("i,j") = Xab("X,i,j") * (V_inv_oh("X,P") * trans("P"));
    }
    J.get_world().gop.fence();

    auto EriK
          = BlockSparseIntegrals(world, eng_pool,
                                 utility::make_array(df_bs, obs, min_bs), op);
    decltype(EriK) Perm, Dterm;
    Dterm("X,a,b") = V_inv_oh("X,P") * EriK("P,a,b");
    Perm("X,b,a") = Dterm("X,a,b");
    K("i,j") = Perm("X,k,i") * (Dterm("X,j,a") * D_min("a,k"));
    F("i,j") = H("i,j") + 2 * J("i,j") - K("i,j");
    decltype(F)::wait_for_lazy_cleanup(F.get_world(), 1);

    return F;
}

} // namespace scf
} // namespace integrals
} // namespace tcc

#endif // TCC_INTEGRALS_SCF_SOAD_H
