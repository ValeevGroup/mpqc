#pragma once
#ifndef TCC_INTEGRALS_SCF_SOAD_H
#define TCC_INTEGRALS_SCF_SOAD_H

#include <mpqc/chemistry/qc/basis/basis.h>
#include <mpqc/chemistry/qc/basis/basis_set.h>
#include <mpqc/chemistry/molecule/molecule.h>

#include "../../../../../include/tiledarray.h"
#include "../../../../../common/typedefs.h"
#include "../../../../../utility/make_array.h"
#include "../../../../../utility/array_info.h"

#include "../../../../../common/namespaces.h"

#include <mpqc/chemistry/qc/integrals/task_integrals.h>
#include <mpqc/chemistry/qc/integrals/task_integrals_common.h>
#include <mpqc/chemistry/qc/integrals/direct_task_integrals.h>

#include "../../../../../ta_routines/array_to_eigen.h"

#include <memory>
#include <vector>

namespace mpqc {
namespace scf {

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

template <typename Engs, typename Array, typename Op>
void soad_task(Engs eng_pool, int64_t ord, ShellVec const *obs_row,
               ShellVec const *obs_col, ShellVec const *min_bs,
               const MatrixD *D, Array *F, Op op) {

    auto range = F->trange().make_tile_range(ord);
    const auto f_extent = range.extent();
    const auto lb = range.lobound();
    TA::TensorD tile(range, 0.0);

    auto &d = *D;

    auto &eng = eng_pool->local();
    eng.set_precision(1e-12);

    using Rng = std::pair<int64_t, int64_t>;

    auto J = [&](double const *buff, Rng sh0_rng, Rng sh1_rng, Rng sh_min_rng) {
        const auto row_start = sh0_rng.first;
        const auto col_start = sh1_rng.first;
        const auto min_start = sh_min_rng.first;

        const auto row_end = sh0_rng.second;
        const auto col_end = sh1_rng.second;
        const auto min_end = sh_min_rng.second;

        const auto row_size = row_end - row_start;
        const auto col_size = col_end - col_start;
        const auto min_size = min_end - min_start;

        for (auto p = 0, j_ord = 0; p < row_size; ++p) {
            const auto fp = p + row_start + lb[0];

            for (auto q = 0; q < col_size; ++q) {
                const auto fq = q + col_start + lb[1];

                auto val = 0.0;
                for (auto r = 0; r < min_size; ++r) {
                    const auto dr = r + min_start;

                    for (auto s = 0; s < min_size; ++s, ++j_ord) {
                        const auto ds = s + min_start;

                        val += buff[j_ord] * d(dr, ds);
                    }
                }

                tile(fp, fq) += 2 * val;
            }
        }
    };

    auto K = [&](double const *buff, Rng sh0_rng, Rng sh1_rng, Rng sh_min_rng) {
        const auto row_start = sh0_rng.first;
        const auto col_start = sh1_rng.first;
        const auto min_start = sh_min_rng.first;

        const auto row_end = sh0_rng.second;
        const auto col_end = sh1_rng.second;
        const auto min_end = sh_min_rng.second;

        const auto row_size = row_end - row_start;
        const auto col_size = col_end - col_start;
        const auto min_size = min_end - min_start;

        for (auto p = 0, k_ord = 0; p < row_size; ++p) {
            const auto fp = p + row_start + lb[0];

            for (auto r = 0; r < min_size; ++r) {
                const auto dr = r + min_start;

                for (auto q = 0; q < col_size; ++q) {
                    const auto fq = q + col_start + lb[1];

                    for (auto s = 0; s < min_size; ++s, ++k_ord) {
                        const auto ds = s + min_start;

                        tile(fp, fq) -= buff[k_ord] * d(dr, ds);
                    }
                }
            }
        }

    };

    // Loop over shells
    auto sh0_start = 0;
    for (auto const &sh0 : *obs_row) {
        const auto sh0_size = sh0.size();
        auto sh0_rng = std::make_pair(sh0_start, sh0_start + sh0_size);
        sh0_start += sh0_size;

        auto sh1_start = 0;
        for (auto const &sh1 : *obs_col) {
            const auto sh1_size = sh1.size();
            auto sh1_rng = std::make_pair(sh1_start, sh1_start + sh1_size);
            sh1_start += sh1_size;

            auto min_start = 0;
            for (auto const &sh_min : *min_bs) {
                const auto min_size = sh_min.size();
                auto min_rng = std::make_pair(min_start, min_start + min_size);
                min_start += min_size;

                const auto& J_bufs = eng.compute(sh0, sh1, sh_min, sh_min);
                TA_USER_ASSERT(J_bufs.size() == 1, "unexpected result from Engine::compute()");
                J(J_bufs[0], sh0_rng, sh1_rng, min_rng);

                const auto& K_bufs = eng.compute(sh0, sh_min, sh1, sh_min);
                TA_USER_ASSERT(K_bufs.size() == 1, "unexpected result from Engine::compute()");
                K(K_bufs[0], sh0_rng, sh1_rng, min_rng);
            }
        }
    }

    F->set(ord, op(std::move(tile)));
    eng.set_precision(integrals::detail::integral_engine_precision);
}

template <typename ShrPool, typename Array,
          typename Op = integrals::TensorPassThrough>
Array fock_from_soad(madness::World &world,
                     molecule::Molecule const &clustered_mol,
                     basis::Basis const &obs, ShrPool engs, Array const &H,
                     Op op = Op{}) {

    // Soad Density
    auto D = soad_density_eig_matrix(clustered_mol);

    // Get minimal basis TODO fix this to BCAST basis to avoid file reads
    const auto min_bs_shells
          = basis::BasisSet("sto-3g").get_flat_shells(clustered_mol);

    // Make F scaffolding
    auto const &trange = H.trange();
    auto const &shape_range = H.get_shape().data().range();

    const auto max_norm = std::numeric_limits<float>::max();
    auto shape_norms = TA::Tensor<float>(shape_range, max_norm);
    TA::SparseShape<float> F_shape(shape_norms, trange);

    Array F(world, trange, F_shape);

    // Loop over lower diagonal tiles
    const auto F_extent = F.trange().tiles().extent();
    for (auto i = 0; i < F_extent[0]; ++i) {
        const auto i_ord = i * F_extent[1];

        for (auto j = 0; j < F_extent[1]; ++j) {
            const auto ord = i_ord + j;

            if (!F.is_zero(ord) && F.is_local(ord)) {
                auto const &obs_row = obs.cluster_shells()[i];
                auto const &obs_col = obs.cluster_shells()[j];

                world.taskq.add(soad_task<ShrPool, Array, Op>, engs, ord,
                                &obs_row, &obs_col, &min_bs_shells, &D, &F, op);
            }
        }
    }
    world.gop.fence();

    F("i,j") += H("i,j");
    F.truncate();
    return F;
}

} // namespace scf
} // namespace mpqc

#endif // TCC_INTEGRALS_SCF_SOAD_H
