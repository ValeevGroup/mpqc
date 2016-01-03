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
    {
        J("i,j") = eri_j("i,j,k,l") * D_min("k,l");
        world.gop.fence();
    }
    {
        K("i,j") = eri_k("i,k,j,l") * D_min("k,l");
        world.gop.fence();
    }
    F("i,j") = H("i,j") + 2 * J("i,j") - K("i,j");
    world.gop.fence();

    return F;
}

} // namespace scf
} // namespace mpqc

#endif // TCC_INTEGRALS_SCF_SOAD_H
