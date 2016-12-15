
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_SOAD_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_SOAD_H_

#include <libint2/chemistry/sto3g_atomic_density.h>

#include "mpqc/chemistry/molecule/molecule.h"
#include "mpqc/chemistry/qc/basis/basis.h"
#include "mpqc/chemistry/qc/basis/basis_set.h"

#include "mpqc/math/external/tiledarray/array_info.h"
#include "mpqc/util/meta/make_array.h"
#include <tiledarray.h>

#include "mpqc/chemistry/qc/integrals/direct_task_integrals.h"
#include "mpqc/chemistry/qc/integrals/task_integrals.h"
#include "mpqc/chemistry/qc/integrals/task_integrals_common.h"

#include "mpqc/math/external/eigen/eigen.h"

#include <memory>
#include <vector>

namespace mpqc {
namespace lcao {

inline RowMatrixXd soad_density_eig_matrix(Molecule const &mol) {
  auto nao = 0;
  for (const auto &atom : mol.atoms()) {
    const auto Z = atom.charge();
    nao += libint2::sto3g_num_ao(Z);
  }

  RowMatrixXd D(nao, nao);
  D.setZero();

  size_t ao = 0;
  for (const auto &atom : mol.atoms()) {
    const auto Z = atom.charge();
    const auto &occvec = libint2::sto3g_ao_occupation_vector(Z);
    for (const auto &occ : occvec) {
      D(ao, ao) = occ;
      ++ao;
    }
  }

  return D * 0.5;  // we use densities normalized to # of electrons/2
}

namespace gaussian {
template <typename Engs, typename Tile, typename Policy>
void soad_task(Engs eng_pool, int64_t ord,
               std::vector<libint2::Shell> const *obs_row,
               std::vector<libint2::Shell> const *obs_col,
               std::vector<libint2::Shell> const *min_bs, const RowMatrixXd *D,
               TA::DistArray<Tile,Policy> *F, std::function<Tile(TA::TensorD &&)> op) {
  auto range = F->trange().make_tile_range(ord);
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

        const auto &J_bufs = eng.compute(sh0, sh1, sh_min, sh_min);
        TA_USER_ASSERT(J_bufs.size() == 1,
                       "unexpected result from Engine::compute()");
        if (J_bufs[0] != nullptr) {
          J(J_bufs[0], sh0_rng, sh1_rng, min_rng);
        }

        const auto &K_bufs = eng.compute(sh0, sh_min, sh1, sh_min);
        TA_USER_ASSERT(K_bufs.size() == 1,
                       "unexpected result from Engine::compute()");
        if (K_bufs[0] != nullptr) {
          K(K_bufs[0], sh0_rng, sh1_rng, min_rng);
        }
      }
    }
  }

  F->set(ord, op(std::move(tile)));
  eng.set_precision(detail::integral_engine_precision);
}

/**
 * fock matrix computed from soad for SparsePolicy
 * @tparam ShrPool
 * @tparam Tile Tile type
 * @tparam Policy TA::SparsePolicy
 * @param world  world object
 * @param clustered_mol  molecule class
 * @param obs basis object
 * @param engs engine
 * @param H  H matrix
 * @param op  operator to convert TA::TensorD to Tile
 * @return Fock matrix from soad
 */
template <typename ShrPool, typename Tile, typename Policy>
TA::DistArray<Tile,typename std::enable_if<std::is_same<Policy, TA::SparsePolicy>::value, TA::SparsePolicy>::type>
fock_from_soad(
    madness::World &world, Molecule const &clustered_mol,
    Basis const &obs, ShrPool engs, TA::DistArray<Tile,Policy> const &H,
    std::function<Tile(TA::TensorD &&)> op = TA::Noop<TA::TensorD, true>()) {
  // Soad Density
  auto D = soad_density_eig_matrix(clustered_mol);

  // Get minimal basis
  const auto min_bs_shells =
      parallel_construct_basis(world, gaussian::BasisSet("sto-3g"), clustered_mol)
          .flattened_shells();
  // Make F scaffolding
  auto const &trange = H.trange();
  auto const &shape_range = H.shape().data().range();

  const auto max_norm = std::numeric_limits<float>::max();
  auto shape_norms = TA::Tensor<float>(shape_range, max_norm);
  TA::SparseShape<float> F_shape(shape_norms, trange);

  TA::DistArray<Tile,Policy> F(world, trange, F_shape);

  // Loop over lower diagonal tiles
  const auto F_extent = F.trange().tiles_range().extent();
  for (auto i = 0; i < F_extent[0]; ++i) {
    const auto i_ord = i * F_extent[1];

    for (auto j = 0; j < F_extent[1]; ++j) {
      const auto ord = i_ord + j;

      if (!F.is_zero(ord) && F.is_local(ord)) {
        auto const &obs_row = obs.cluster_shells()[i];
        auto const &obs_col = obs.cluster_shells()[j];

        world.taskq.add(soad_task<ShrPool,Tile,Policy>, engs, ord, &obs_row,
                        &obs_col, &min_bs_shells, &D, &F, op);
      }
    }
  }
  world.gop.fence();

  F("i,j") += H("i,j");
  F.truncate();
  return F;
}


/**
 * fock matrix computed from soad for DensePolicy
 * @tparam ShrPool
 * @tparam Tile Tile type
 * @tparam Policy TA::DensePolicy
 * @param world  world object
 * @param clustered_mol  molecule class
 * @param obs basis object
 * @param engs engine
 * @param H  H matrix
 * @param op  operator to convert TA::TensorD to Tile
 * @return Fock matrix from soad
 */
template <typename ShrPool, typename Tile, typename Policy>
TA::DistArray<Tile,typename std::enable_if<std::is_same<Policy, TA::DensePolicy>::value, TA::DensePolicy>::type>
fock_from_soad(
    madness::World &world, Molecule const &clustered_mol,
    Basis const &obs, ShrPool engs, TA::DistArray<Tile,Policy> const &H,
    std::function<Tile(TA::TensorD &&)> op = TA::Noop<TA::TensorD, true>()) {
  // Soad Density
  auto D = soad_density_eig_matrix(clustered_mol);

  // Get minimal basis
  const auto min_bs_shells =
      parallel_construct_basis(world, gaussian::BasisSet("sto-3g"), clustered_mol)
          .flattened_shells();

  auto const &trange = H.trange();
  TA::DistArray<Tile,Policy> F(world, trange);

  // Loop over lower diagonal tiles
  const auto F_extent = F.trange().tiles_range().extent();
  for (auto i = 0; i < F_extent[0]; ++i) {
    const auto i_ord = i * F_extent[1];

    for (auto j = 0; j < F_extent[1]; ++j) {
      const auto ord = i_ord + j;

      if (F.is_local(ord)) {
        auto const &obs_row = obs.cluster_shells()[i];
        auto const &obs_col = obs.cluster_shells()[j];

        world.taskq.add(soad_task<ShrPool,Tile, Policy>, engs, ord, &obs_row,
                        &obs_col, &min_bs_shells, &D, &F, op);
      }
    }
  }
  world.gop.fence();

  F("i,j") += H("i,j");
  return F;
}


}  // namespace gaussian
}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_SOAD_H_
