#pragma once
#ifndef MPQC_SCF_CADFFITTINGCOEFFS_H
#define MPQC_SCF_CADFFITTINGCOEFFS_H

#include "../include/tiledarray.h"
#include "../common/namespaces.h"
#include "../integrals/integrals.h"
#include "../common/typedefs.h"

#include "../molecule/molecule.h"
#include "../basis/basis.h"
#include "../basis/basis_set.h"

#include "../utility/make_array.h"

#include <unordered_map>

namespace mpqc {
namespace scf {

namespace cadf {

inline basis::Basis by_atom_basis(
      molecule::Molecule const &mol, basis::BasisSet const &bs,
      std::unordered_map<std::size_t, std::size_t> &atom_to_cluster_map) {

    std::vector<molecule::AtomBasedClusterable> atoms;
    std::unordered_map<std::size_t, std::size_t> obs_atom_to_cluster_map;

    auto cluster_ind = 0;
    auto atom_ind = 0;
    for (auto const &cluster : mol.clusterables()) {
        for (auto atom : cluster.atoms()) {
            atoms.push_back(std::move(atom));
            atom_to_cluster_map[atom_ind] = cluster_ind;
            ++atom_ind;
        }
        ++cluster_ind;
    }

    auto sort_atoms_in_mol = false;
    return basis::Basis(bs.get_cluster_shells(
          molecule::Molecule(std::move(atoms), sort_atoms_in_mol)));
}


inline TA::TiledRange
cadf_trange(basis::Basis const &obs_by_atom, basis::Basis const &dfbs_by_atom) {
    std::vector<TA::TiledRange1> trange1s;
    trange1s.emplace_back(dfbs_by_atom.create_trange1());
    trange1s.emplace_back(obs_by_atom.create_trange1());
    trange1s.emplace_back(obs_by_atom.create_trange1());

    return TRange(trange1s.begin(), trange1s.end());
}

inline TA::SparseShape<float>
cadf_shape(madness::World &world, basis::Basis const &obs,
           TA::TiledRange const &trange) {

    auto &tiles = trange.tiles();
    TA::TensorD norms(tiles, 0.0);

    for (auto i = 0ul; i < obs.nclusters(); ++i) {
        for (auto j = 0ul; j < obs.nclusters(); ++j) {
            auto idx0 = std::array<unsigned long, 3>{{i, i, j}};
            auto idx1 = std::array<unsigned long, 3>{{j, i, j}};

            auto ord0 = tiles.ordinal(idx0);
            auto ord1 = tiles.ordinal(idx1);

            *(norms.data() + ord0) = std::numeric_limits<double>::max();
            *(norms.data() + ord1) = std::numeric_limits<double>::max();
        }
    }

    return TA::SparseShape<float>(world, norms, trange);
}

template <typename Integral>
void create_tiles(madness::World &world,
                  TA::DistArray<TA::TensorD, SpPolicy> &C_df, Integral &eri3,
                  TA::DistArray<TA::TensorD, TA::DensePolicy> const &M,
                  unsigned long natoms) {

    auto task = [&](unsigned long i, unsigned long j) {
        auto &trange = C_df.trange();
        auto &tiles = trange.tiles();
        auto &M_tiles = M.trange().tiles();

        if (i == j) {
            auto eri3_idx = std::array<unsigned long, 3>{{i, i, i}};
            auto eri3_ord = tiles.ordinal(eri3_idx);

            if (C_df.is_local(eri3_ord)) {
                TA::TensorD eri3_tile = eri3.array().find(eri3_ord).get();
                auto eri3_extent = eri3_tile.range().extent();
                MatrixD eri3_eig
                      = TA::eigen_map(eri3_tile, eri3_extent[0],
                                      eri3_extent[1] * eri3_extent[2]);

                auto M_idx = std::array<unsigned long, 2>{{i, i}};
                auto M_ord = M_tiles.ordinal(M_idx);
                auto M_tile = M.find(M_ord).get();
                auto M_extent = M_tile.range().extent();
                MatrixD M_eig = TA::eigen_map(M_tile, M_extent[0], M_extent[1]);

                TA::TensorD out_tile(eri3_tile.range());

                MatrixD out_eig = M_eig.inverse() * eri3_eig;
                TA::eigen_map(out_tile, eri3_extent[0],
                              eri3_extent[1] * eri3_extent[2]) = out_eig;

                C_df.set(eri3_ord, out_tile);
            }
        } else {
            auto eri3_idx_0 = std::array<unsigned long, 3>{{i, i, j}};
            auto eri3_ord_0 = tiles.ordinal(eri3_idx_0);

            auto eri3_idx_1 = std::array<unsigned long, 3>{{j, i, j}};
            auto eri3_ord_1 = tiles.ordinal(eri3_idx_1);

            if (C_df.is_local(eri3_idx_0) || C_df.is_local(eri3_idx_1)) {

                auto eri3_range_0 = trange.make_tile_range(eri3_idx_0);
                auto eri3_range_1 = trange.make_tile_range(eri3_idx_1);

                auto eri3_extent_0 = eri3_range_0.extent();
                auto eri3_extent_1 = eri3_range_1.extent();

                TA::TensorD eri3_tile_0 = eri3.array().find(eri3_ord_0).get();
                MatrixD eri3_eig_0
                      = TA::eigen_map(eri3_tile_0, eri3_extent_0[0],
                                      eri3_extent_0[1] * eri3_extent_0[2]);

                TA::TensorD eri3_tile_1 = eri3.array().find(eri3_ord_1).get();
                MatrixD eri3_eig_1
                      = TA::eigen_map(eri3_tile_1, eri3_extent_1[0],
                                      eri3_extent_1[1] * eri3_extent_1[2]);

                auto M_idx_11 = std::array<unsigned long, 2>{{i, i}};
                auto M_idx_12 = std::array<unsigned long, 2>{{i, j}};
                auto M_idx_21 = std::array<unsigned long, 2>{{j, i}};
                auto M_idx_22 = std::array<unsigned long, 2>{{j, j}};

                auto M_ord_11 = M_tiles.ordinal(M_idx_11);
                auto M_ord_12 = M_tiles.ordinal(M_idx_12);
                auto M_ord_21 = M_tiles.ordinal(M_idx_21);
                auto M_ord_22 = M_tiles.ordinal(M_idx_22);

                TA::TensorD M_tile_11 = M.find(M_ord_11).get();
                TA::TensorD M_tile_12 = M.find(M_ord_12).get();
                TA::TensorD M_tile_21 = M.find(M_ord_21).get();
                TA::TensorD M_tile_22 = M.find(M_ord_22).get();

                auto M_extent_11 = M_tile_11.range().extent();
                auto M_extent_12 = M_tile_12.range().extent();
                auto M_extent_21 = M_tile_21.range().extent();
                auto M_extent_22 = M_tile_22.range().extent();

                auto rows = M_extent_11[0] + M_extent_21[0];
                auto cols = M_extent_11[1] + M_extent_12[1];

                MatrixD M_combo(rows, cols);
                // Write M11
                M_combo.block(0, 0, M_extent_11[0], M_extent_11[1])
                      = TA::eigen_map(M_tile_11, M_extent_11[0],
                                      M_extent_11[1]);

                // Write M12
                M_combo.block(0, M_extent_11[1], M_extent_12[0], M_extent_12[1])
                      = TA::eigen_map(M_tile_12, M_extent_12[0],
                                      M_extent_12[1]);

                // Write M21
                M_combo.block(M_extent_11[0], 0, M_extent_21[0], M_extent_21[1])
                      = TA::eigen_map(M_tile_21, M_extent_21[0],
                                      M_extent_21[1]);

                // Write M22
                M_combo.block(M_extent_11[0], M_extent_11[1], M_extent_22[0],
                              M_extent_22[1])
                      = TA::eigen_map(M_tile_22, M_extent_22[0],
                                      M_extent_22[1]);

                MatrixD M_combo_inv = M_combo.inverse();

                // Doing the block wise GEMM by hand for now.
                auto block11
                      = M_combo_inv.block(0, 0, M_extent_11[0], M_extent_11[1]);
                auto block12
                      = M_combo_inv.block(0, M_extent_11[1], M_extent_12[0],
                                          M_extent_12[1]);
                auto block21
                      = M_combo_inv.block(M_extent_11[0], 0, M_extent_21[0],
                                          M_extent_21[1]);
                auto block22
                      = M_combo_inv.block(M_extent_11[0], M_extent_11[1],
                                          M_extent_22[0], M_extent_22[1]);

                MatrixD C0 = block11 * eri3_eig_0 + block12 * eri3_eig_1;
                MatrixD C1 = block21 * eri3_eig_0 + block22 * eri3_eig_1;

                TA::TensorD out_tile_0(eri3_range_0);
                TA::TensorD out_tile_1(eri3_range_1);

                TA::eigen_map(out_tile_0, eri3_extent_0[0],
                              eri3_extent_0[1] * eri3_extent_0[2]) = C0;

                TA::eigen_map(out_tile_1, eri3_extent_1[0],
                              eri3_extent_1[1] * eri3_extent_1[2]) = C1;

                if (C_df.is_local(eri3_ord_0)) {
                    C_df.set(eri3_ord_0, out_tile_0);
                }

                if (C_df.is_local(eri3_ord_1)) {
                    C_df.set(eri3_ord_1, out_tile_1);
                }
            }
        }
    };

    for (auto i = 0ul; i < natoms; ++i) {
        for (auto j = 0ul; j < natoms; ++j) {
            world.taskq.add(task, i, j);
        }
    }
    world.gop.fence();
}

} // namespace cadf


inline TA::DistArray<TA::TensorD, SpPolicy> compute_atomic_fitting_coeffs(
      madness::World &world, molecule::Molecule const &obs_molecule,
      molecule::Molecule const &dfbs_molecule, basis::BasisSet const &obs_set,
      basis::BasisSet const &dfbs_set,
      std::unordered_map<std::size_t, std::size_t> &obs_atom_to_cluster_map,
      std::unordered_map<std::size_t, std::size_t> &dfbs_atom_to_cluster_map) {

    auto by_atom_obs
          = cadf::by_atom_basis(obs_molecule, obs_set, obs_atom_to_cluster_map);

    auto by_atom_dfbs = cadf::by_atom_basis(dfbs_molecule, dfbs_set,
                                            dfbs_atom_to_cluster_map);

    const auto dfbs_array = utility::make_array(by_atom_dfbs, by_atom_dfbs);
    auto eri_e = integrals::make_2body_shr_pool(by_atom_dfbs, by_atom_obs);

    auto M = integrals::dense_integrals(world, eri_e, dfbs_array);

    auto three_c_array
          = utility::make_array(by_atom_dfbs, by_atom_obs, by_atom_obs);

    auto eri3 = integrals::untruncated_direct_sparse_integrals(world, eri_e,
                                                               three_c_array);

    auto trange = cadf::cadf_trange(by_atom_obs, by_atom_dfbs);

    auto shape = cadf::cadf_shape(world, by_atom_obs, trange);

    auto pmap = eri3.array().get_pmap();

    TA::DistArray<TA::TensorD, SpPolicy> C_df(world, trange, shape, pmap);

    cadf::create_tiles(world, C_df, eri3, M, by_atom_obs.nclusters());
    C_df.truncate();
    world.gop.fence();

    return C_df;
}

} // namespace scf
} // namespace mpqc

#endif // MPQC_SCF_CADFFITTINGCOEFFS_H
