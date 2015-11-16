#pragma once
#ifndef TCC_INTEGRALS_SCF_SOAD_H
#define TCC_INTEGRALS_SCF_SOAD_H

#include "../basis/basis.h"
#include "../basis/basis_set.h"
#include "../molecule/molecule.h"

#include "../include/tiledarray.h"
#include "../common/typedefs.h"
#include "../utility/make_array.h"
#include "../utility/array_storage.h"

#include "../common/namespaces.h"

#include "../integrals/task_integrals.h"
#include "../integrals/task_integrals_common.h"
#include "../integrals/direct_task_integrals.h"

#include "../ta_routines/array_to_eigen.h"

#include <memory>
#include <vector>

namespace mpqc {
namespace scf {

template <typename Op>
using soad_tile = integrals::detail::Ttype<Op>;

/*! \brief constructs a tile using Edwards SOAD strategy.
 *
 * This function creates a tile to initialize a density guess based on
 * the Superposition of Atomic Densities (SOAD) approach.  The function
 * expects to receive symmetric diagonal ranges.
 */
template <typename Op = integrals::TensorPassThrough>
soad_tile<Op>
soad_tile_build(TA::Range range, molecule::AtomBasedClusterable cluster,
                Op op = Op{}) {

    TA::Tensor<double> tensor(std::move(range), 0.0);
    auto const &extent = tensor.range().extent();
    auto t_map = TA::eigen_map(tensor, extent[0], extent[1]);

    auto atoms = cluster.atoms();
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
            const double num_electrons_per_2p
                  = (Z > 4) ? (double)(Z - 4) / 3 : 0;
            for (auto xyz = 0; xyz != 3; ++xyz)
                t_map(ao_offset + 2 + xyz, ao_offset + 2 + xyz)
                      = num_electrons_per_2p;
            ao_offset += 5;
        }
    }

    // Multiply by 1/2 because we are using the convention that the Density is
    // idempotent.
    t_map *= 0.5;

    return op(std::move(tensor));
}

template <typename Op = integrals::TensorPassThrough>
DArray<2, soad_tile<Op>, SpPolicy>
soad_density_guess(madness::World &world, molecule::Molecule clustered_mol,
                   basis::Basis const &min_bs, Op op = Op{}) {

    auto trange1 = min_bs.create_trange1();
    auto trange = TRange({trange1, trange1});

    TA::Tensor<float> tile_norms(trange.tiles(), 0.0);

    // Unzero diagonal tiles may break for non-symetric tiling
    auto const extent = tile_norms.range().extent();
    auto tn_map = TA::eigen_map(tile_norms, extent[0], extent[1]);
    for (auto i = 0; i < tn_map.rows(); ++i) {
        tn_map(i, i) = std::numeric_limits<float>::max(); // Go big
    }

    // Create the array and shape
    TA::SparseShape<float> shape(world, tile_norms, trange);
    DArray<2, soad_tile<Op>, SpPolicy> D_min(world, trange, shape);

    // Initialize D_min tiles by calling soad_tile
    for (auto it : *D_min.get_pmap()) {
        if (!D_min.is_zero(it) && D_min.is_local(it)) {
            auto idx = trange.tiles().idx(it);
            auto range = trange.make_tile_range(it);
            if (idx[0] == idx[1]) {
                madness::Future<soad_tile<Op>> tile
                      = world.taskq.add(soad_tile_build<Op>, range,
                                        clustered_mol.clusterables()[idx[0]],
                                        op);
                D_min.set(it, tile);
            } else {

                auto empty = [](TA::Range rng, Op op) {
                    return op(TA::TensorD(rng, 0.0));
                };

                madness::Future<soad_tile<Op>> tile
                      = world.taskq.add(empty, range, op);
                D_min.set(it, tile);
            }
        }
    }
    world.gop.fence();

    D_min.truncate(); // Fix tile norms.
    return D_min;
}

template <typename ShrPool, typename Array,
          typename Op = integrals::TensorPassThrough>
Array fock_from_soad(madness::World &world, molecule::Molecule clustered_mol,
                     basis::Basis const &obs, ShrPool engs, Array const &H,
                     Op op = Op{}) {

    auto min_bs = basis::Basis(
          basis::BasisSet("sto-3g").get_cluster_shells(clustered_mol));

    auto obs_trange1 = obs.create_trange1();
    auto minbs_trange1 = min_bs.create_trange1();

    auto trange_j
          = TRange({obs_trange1, obs_trange1, minbs_trange1, minbs_trange1});
    auto trange_k
          = TRange({obs_trange1, minbs_trange1, obs_trange1, minbs_trange1});

    auto bs_j = tcc::utility::make_array(obs, obs, min_bs, min_bs);
    auto bs_k = tcc::utility::make_array(obs, min_bs, obs, min_bs);

    auto eri_j = integrals::soad_direct_integrals(world, engs, bs_j);
    auto eri_k = integrals::soad_direct_integrals(world, engs, bs_k);

    auto D_min = soad_density_guess(world, clustered_mol, min_bs, op);

    Array F, J, K;
    J("i,j") = eri_j("i,j,k,l") * D_min("k,l");
    K("i,j") = eri_k("i,k,j,l") * D_min("k,l");
    F("i,j") = H("i,j") + 2 * J("i,j") - K("i,j");

    return F;
}

} // namespace scf
} // namespace mpqc

#if 0

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
    // make range for decomposed tensor
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
            const double num_electrons_per_2p
                  = (Z > 4) ? (double)(Z - 4) / 3 : 0;
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

/*! \brief just construsts a tile of zero's
*/
TileType empty_tile(TA::Range range, double cut) {
    // make range for decomposed tensor
    const auto i = range.extent()[0];
    const auto j = range.extent()[1];

    TA::Range local_range{i, j};
    TA::Tensor<double> tensor(std::move(local_range), 0.0);

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
            auto idx = trange.tiles().idx(*it);
            if (idx[0] == idx[1]) {
                auto cluster_ord = trange.tiles().idx(*it)[0];
                madness::Future<tensor::Tile<tensor::DecomposedTensor<double>>>
                      tile = world.taskq.add(soad_tile, clusters[cluster_ord],
                                             trange.make_tile_range(*it), cut);
                D_min.set(*it, tile);
            } else {
                madness::Future<tensor::Tile<tensor::DecomposedTensor<double>>>
                      tile = world.taskq.add(empty_tile,
                                             trange.make_tile_range(*it), cut);
                D_min.set(*it, tile);
            }
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

    decltype(D_min) L_d;
    {
        auto Dm_eig = array_ops::array_to_eigen(D_min);
        Eig::LLT<decltype(Dm_eig)> llt(Dm_eig);
        decltype(Dm_eig) L_d_eig = llt.matrixL();
        L_d = array_ops::eigen_to_array<decltype(D_min)::value_type>(
              D_min.get_world(), L_d_eig, D_min.trange().data()[0],
              D_min.trange().data()[1], cut);
    }


    decltype(D_min) J, K, F;

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
    // Reuse EriK so we don't have a temp Array laying around.
    EriK("X,i,a") = V_inv_oh("X,P") * (EriK("P,a,b") * L_d("b,i"));
    K("i,j") = EriK("X,k,i") * EriK("X,k,j");
    F("i,j") = H("i,j") + 2 * J("i,j") - K("i,j");
    F.get_world().gop.fence();

    return F;
}

template <typename SharedEnginePool, typename Op>
ArrayType fock_from_minimal_low_mem(
      madness::World &world, basis::Basis const &obs, basis::Basis const &df_bs,
      SharedEnginePool eng_pool, ArrayType const &H, ArrayType const &V_inv_oh,
      ArrayType const &V_inv, Eri3ArrayType const &Xab,
      std::vector<std::shared_ptr<molecule::Cluster>> const &clusters,
      double cut, Op op) {

    basis::BasisSet min_bs_set("sto-3g");

    std::streambuf *cout_sbuf = std::cout.rdbuf(); // Silence libint printing.
    std::ofstream fout("/dev/null");
    std::cout.rdbuf(fout.rdbuf());
    basis::Basis min_bs{min_bs_set.create_soad_basis(clusters)};
    std::cout.rdbuf(cout_sbuf);
    auto D_min = minimal_density_guess(world, clusters, min_bs, cut);

    decltype(D_min) L_d;
    {
        auto Dm_eig = array_ops::array_to_eigen(D_min);
        Eig::LLT<decltype(Dm_eig)> llt(Dm_eig);
        decltype(Dm_eig) L_d_eig = llt.matrixL();
        L_d = array_ops::eigen_to_array<decltype(D_min)::value_type>(
              D_min.get_world(), L_d_eig, D_min.trange().data()[0],
              D_min.trange().data()[1], cut);
    }


    decltype(D_min) J, K, F;
    {
        auto EriJ = BlockSparseIntegrals(
              world, eng_pool, utility::make_array(df_bs, min_bs, min_bs), op);
        TA::Array<double, 1, typename decltype(EriJ)::value_type,
                  TA::SparsePolicy> trans;
        trans("P") = V_inv("P,X") * EriJ("X,a,b") * D_min("a,b");
        J("i,j") = Xab("X,i,j") * (trans("P"));
    }
    J.get_world().gop.fence();

    auto EriK
          = BlockSparseIntegrals(world, eng_pool,
                                 utility::make_array(df_bs, obs, min_bs), op);
    // Reuse EriK so we don't have a temp Array laying around.
    EriK("X,i,a") = V_inv_oh("X,P") * (EriK("P,a,b") * L_d("b,i"));
    K("i,j") = EriK("X,k,i") * EriK("X,k,j");
    F("i,j") = H("i,j") + 2 * J("i,j") - K("i,j");
    F.get_world().gop.fence();

    return F;
}

} // namespace scf
} // namespace integrals
} // namespace tcc
#endif

#endif // TCC_INTEGRALS_SCF_SOAD_H
