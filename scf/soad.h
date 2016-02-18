#pragma once
#ifndef TCC_INTEGRALS_SCF_SOAD_H
#define TCC_INTEGRALS_SCF_SOAD_H

#include "../basis/basis.h"
#include "../basis/basis_set.h"
#include "../molecule/molecule.h"

#include "../include/tiledarray.h"
#include "../common/typedefs.h"
#include "../utility/make_array.h"
#include "../utility/array_info.h"

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

    return op(std::move(tensor));
}

template <typename Op = integrals::TensorPassThrough>
DArray<2, soad_tile<Op>, SpPolicy>
soad_density_guess(madness::World &world, molecule::Molecule clustered_mol,
                   basis::Basis const &min_bs, Op op = Op{}) {

    auto trange1 = min_bs.create_trange1();
    auto trange = TRange({trange1, trange1});

    TA::Tensor<float> tile_norms(trange.tiles(), 0.0);

    // non-zero diagonal tiles may break for non-symetric tiling
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

MatrixD soad_density_eig_matrix(molecule::Molecule const &mol) {

    auto nao = 0;
    for (const auto &atom : mol.atoms()) {
        const auto Z = atom.charge();
        if (Z == 1 || Z == 2) // H, He
            nao += 1;
        else if (Z <= 10) // Li - Ne
            nao += 5;
        else
            throw "SOAD with Z > 10 is not yet supported";
    }

    MatrixD D(nao, nao);
    D.setZero();

    size_t ao_offset = 0; // first AO of this atom
    for (const auto &atom : mol.atoms()) {
        const auto Z = atom.charge();
        if (Z == 1 || Z == 2) {          // H, He
            D(ao_offset, ao_offset) = Z; // all electrons go to the 1s
            ao_offset += 1;
        } else if (Z <= 10) {
            D(ao_offset, ao_offset) = 2; // 2 electrons go to the 1s
            D(ao_offset + 1, ao_offset + 1)
                  = (Z == 3) ? 1
                             : 2; // Li? only 1 electron in 2s, else 2 electrons
            // smear the remaining electrons in 2p orbitals
            const double num_electrons_per_2p = (Z > 4) ? (double)(Z - 4) / 3
                                                        : 0;
            for (auto xyz = 0; xyz != 3; ++xyz)
                D(ao_offset + 2 + xyz, ao_offset + 2 + xyz)
                      = num_electrons_per_2p;
            ao_offset += 5;
        }
    }

    return D * 0.5; // we use densities normalized to # of electrons/2
}

template <typename ShrPool, typename Array,
          typename Op = integrals::TensorPassThrough>
Array fock_from_soad(madness::World &world,
                     molecule::Molecule const &clustered_mol,
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

    auto bs_j = utility::make_array(obs, obs, min_bs, min_bs);
    auto bs_k = utility::make_array(obs, min_bs, obs, min_bs);

    world.gop.fence();
    auto old_thresh = TA::SparseShape<float>::threshold();

    const auto soad_thresh = 1e-7;
    TiledArray::SparseShape<float>::threshold(soad_thresh);

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
    TiledArray::SparseShape<float>::threshold(old_thresh);
    world.gop.fence();

    return F;
}

template <typename ShrPool, typename Array,
          typename Op = integrals::TensorPassThrough>
Array fock_from_soad_low_mem(madness::World &world,
                             molecule::Molecule const &clustered_mol,
                             basis::Basis const &obs, ShrPool engs,
                             Array const &H, Op op = Op{}) {

    auto D = soad_density_eig_matrix(clustered_mol);

    MatrixD G(obs.nfunctions(), obs.nfunctions());
    G.setZero();

    auto min_bs = basis::Basis(
          basis::BasisSet("sto-3g").get_cluster_shells(clustered_mol));

    auto shell_2func = [](std::vector<libint2::Shell> const &shells) {
        std::vector<int64_t> sh_2_f;

        auto start = 0;
        for (auto const &shell : shells) {
            sh_2_f.push_back(start);
            start += shell.size();
        }

        return sh_2_f;
    };

    const auto nshells = obs.nshells();
    const auto min_nshells = min_bs.nshells();

    const auto obs_shells = obs.flattened_shells();
    const auto shell2bf = shell_2func(obs_shells);

    const auto min_shells = min_bs.flattened_shells();
    const auto shell2bf_D = shell_2func(min_shells);

    auto lambda = [&](decltype(engs) &eng) {
        auto &engine = eng->local();
        engine.set_precision(0);
        auto &g = G;

        for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
            auto const &sh1 = obs_shells[s1];
            auto bf1_first = shell2bf[s1];
            auto n1 = sh1.size();

            for (auto s2 = 0; s2 <= s1; ++s2) {
                auto const &sh2 = obs_shells[s2];
                auto bf2_first = shell2bf[s2];
                auto n2 = sh2.size();

                for (auto s3 = 0; s3 < min_nshells; ++s3) {
                    auto const &sh3 = min_shells[s3];
                    auto bf3_first = shell2bf_D[s3];
                    const auto n3 = sh3.size();

                    for (auto s4 = s3; s4 != s3 + 1; ++s4, ++s1234) {
                        auto const &sh4 = min_shells[s4];
                        auto bf4_first = shell2bf_D[s4];
                        auto n4 = sh4.size();

                        // compute the permutational degeneracy (i.e. # of
                        // equivalents) of the given shell set
                        auto s12_deg = (s1 == s2) ? 1.0 : 2.0;

                        if (s3 >= s4) {
                            auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
                            auto s1234_deg = s12_deg * s34_deg;

                            const auto *buf_J
                                  = engine.compute(sh1, sh2, sh3, sh4);

                            for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
                                const auto bf1 = f1 + bf1_first;
                                for (auto f2 = 0; f2 != n2; ++f2) {
                                    const auto bf2 = f2 + bf2_first;
                                    for (auto f3 = 0; f3 != n3; ++f3) {
                                        const auto bf3 = f3 + bf3_first;
                                        for (auto f4 = 0; f4 != n4;
                                             ++f4, ++f1234) {
                                            const auto bf4 = f4 + bf4_first;

                                            const auto value = buf_J[f1234];
                                            const auto value_scal_by_deg
                                                  = value * s1234_deg;
                                            g(bf1, bf2) += 2.0 * D(bf3, bf4)
                                                           * value_scal_by_deg;
                                        }
                                    }
                                }
                            }
                        }
                        

                        const auto *buf_K = engine.compute(sh1, sh3, sh2, sh4);

                        for (auto f1 = 0, f1324 = 0; f1 != n1; ++f1) {
                            const auto bf1 = f1 + bf1_first;
                            for (auto f3 = 0; f3 != n3; ++f3) {
                                const auto bf3 = f3 + bf3_first;
                                for (auto f2 = 0; f2 != n2; ++f2) {
                                    const auto bf2 = f2 + bf2_first;
                                    for (auto f4 = 0; f4 != n4; ++f4, ++f1324) {
                                        const auto bf4 = f4 + bf4_first;

                                        const auto value = buf_K[f1324];
                                        const auto value_scal_by_deg
                                              = value * s12_deg;
                                        g(bf1, bf2) -= D(bf3, bf4)
                                                       * value_scal_by_deg;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

    }; // thread lambda

    lambda(engs);

    // Have to watch out for aliasing issues
    MatrixD Gt = 0.5 * (G + G.transpose());

    auto obs_trange1 = obs.create_trange1();
    auto G_TA
          = array_ops::eigen_to_array<TA::TensorD>(H.get_world(), Gt,
                                                   obs_trange1, obs_trange1);

    decltype(G_TA) F;
    F("i,j") = H("i,j") + G_TA("i,j");

    world.gop.fence();
    return F;
}

} // namespace scf
} // namespace mpqc

#endif // TCC_INTEGRALS_SCF_SOAD_H
