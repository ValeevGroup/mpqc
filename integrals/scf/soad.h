#pragma once
#ifndef TCC_INTEGRALS_SCF_SOAD_H
#define TCC_INTEGRALS_SCF_SOAD_H

#include "../../basis/basis.h"
#include "../../basis/basis_set.h"
#include "../../molecule/cluster.h"
#include "../../molecule/cluster_collapse.h"

#include "../../include/tiledarray.h"
#include "../../common/typedefs.h"
#include "../../utility/make_array.h"
#include "../../tensor/tcc_tile.h"
#include "../../tensor/decomposed_tensor.h"
#include "../../tensor/decomposed_tensor_nonintrusive_interface.h"

#include "../sparse_task_integrals.h"

#include <memory>
#include <vector>

namespace tcc {
namespace integrals {
namespace scf {

using DataType = tensor::DecomposedTensor<double>;
using TileType = tensor::Tile<DataType>;
using ArrayType = TA::Array<double, 2, TileType, TA::SparsePolicy>;
using Eri3ArrayType = TA::Array<double, 3, TileType, TA::SparsePolicy>;

TileType soad_tile(std::shared_ptr<molecule::Cluster> cluster, TA::Range range,
                   double cut) {
    // make range for decomposed tensor
    const auto i = range.size()[0];
    const auto j = range.size()[1];
    TA::Range local_range{i, j};
    TA::Tensor<double> tensor(std::move(local_range), 0.0);
    auto atoms = molecule::collapse_to_atoms(*cluster);

    auto t_map = TA::eigen_map(tensor, range.size()[0], range.size()[1]);
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
    t_map *= 0.5;
    return TileType(range, DataType(cut, std::move(tensor)));
}


ArrayType minimal_density_guess(
      madness::World &world,
      std::vector<std::shared_ptr<molecule::Cluster>> const &clusters,
      basis::Basis const &min_bs, double cut) {

    TA::TiledRange trange
          = sparse::create_trange(utility::make_array(min_bs, min_bs));
    TA::Tensor<float> tile_norms(trange.tiles(), 0.0);
    auto tn_map = TA::eigen_map(tile_norms, tile_norms.range().size()[0],
                                tile_norms.range().size()[1]);

    auto ord = 0ul;
    for (auto i = 0; i < tn_map.rows(); ++i) {
        auto range = trange.make_tile_range(ord);
        tn_map(i, i) = range.size()[1];
        ord += trange.tiles().weight()[1] + 1;
    }

    TA::SparseShape<float> shape(world, tile_norms, trange);

    ArrayType D_min(world, trange, shape);

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
    basis::Basis min_bs{min_bs_set.create_basis(clusters)};
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
    decltype(D_min) temp;
    temp("i,j") = 2 * J("i,j") - K("i,j");
    F("i,j") = H("i,j") + temp("i,j");
    // F("i,j") = H("i,j") + 2 * J("i,j") - K("i,j");
//    auto tnorm = temp("i,j").norm().get();
    /* auto jnorm = J("i,j").norm().get(); */
    /* auto knorm = K("i,j").norm().get(); */
    /* auto hnorm = H("i,j").norm().get(); */
    /* auto fnorm = F("i,j").norm().get(); */

    /* if (J.get_world().rank() == 0) { */
    /*     std::cout << "Norm J_min = " << jnorm << std::endl; */
    /*     std::cout << "Norm K_min = " << knorm << std::endl; */
    /*     std::cout << "Norm H_min = " << hnorm << std::endl; */
    /*     std::cout << "Norm F_min = " << fnorm << std::endl; */
    /*     std::cout << "Norm temp = " << tnorm << std::endl; */
    /* } */

    return F;
}


} // namespace scf
} // namespace integrals
} // namespace tcc

#endif // TCC_INTEGRALS_SCF_SOAD_H
