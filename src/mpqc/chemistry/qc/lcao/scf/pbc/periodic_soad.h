#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_SOAD_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_SOAD_H_

#include <libint2/chemistry/sto3g_atomic_density.h>
#include <tiledarray.h>

#include "mpqc/chemistry/molecule/unit_cell.h"
#include "mpqc/chemistry/qc/lcao/basis/basis.h"
#include "mpqc/chemistry/qc/lcao/factory/periodic_ao_factory.h"
#include "mpqc/chemistry/qc/lcao/scf/soad.h"

namespace mpqc {
namespace lcao {
namespace gaussian {

/*!
 * \brief This computes SOAD guess for the density matrix in minimal basis
 * for the reference cell, and projects to normal basis of the periodic system
 *
 * \param world MADNESS world
 * \param unitcell UnitCell object
 * \param H core Hamiltonian in real space, should be a complex tile
 * \param ao_factory PeriodicAOFactory object
 * \param op a functor that takes TA::TensorZ && and returns a valid tile type
 * \return Fock matrix in real space
 */
template <typename Tile, typename Policy, typename FactoryType>
TA::DistArray<Tile, Policy> periodic_fock_soad(
    madness::World &world, UnitCell const &unitcell,
    TA::DistArray<Tile, Policy> const &H, FactoryType &pao_factory) {
  ExEnv::out0() << "\nBuilding Fock Matrix from SOAD Guess ...\n";

  using TArray = typename FactoryType::TArray;
  using DirectTArray = typename FactoryType::DirectTArray;

  // get necessary information for periodic AO integrals
  auto RJ_size = pao_factory.RJ_size();
  auto RJ_max = pao_factory.RJ_max();
  auto dcell = unitcell.dcell();
  auto screen = pao_factory.screen();
  auto screen_thresh = pao_factory.screen_threshold();

  // declare necessary variables for periodic AO integrals
  auto g_J_vector = std::vector<DirectTArray>(RJ_size, DirectTArray());
  auto g_K_vector = std::vector<DirectTArray>(RJ_size, DirectTArray());
  std::shared_ptr<utility::TSPool<libint2::Engine>> engine_pool;
  BasisVector bases;
  auto p_screener = std::make_shared<Screener>(Screener{});

  // soad density
  auto D_real = soad_density_eig_matrix(unitcell);
  auto D_comp = D_real.cast<std::complex<double>>();

  // get minimal basis
  auto min_bs = parallel_make_basis(world, Basis::Factory("sto-3g"), unitcell);

  // transform soad density from Eigen to TA
  auto min_bases = BasisVector{{min_bs, min_bs}};
  auto min_trange = detail::create_trange(min_bases);
  auto min_tr0 = min_trange.data()[0];
  auto min_tr1 = min_trange.data()[1];
  auto D =
      array_ops::eigen_to_array<Tile, Policy>(world, D_comp, min_tr0, min_tr1);

  // get normal basis
  Vector3d zero_shift_base(0.0, 0.0, 0.0);
  auto R_max = pao_factory.R_max();
  auto normal_bs = *pao_factory.basis_registry()->retrieve(OrbitalIndex(L"λ"));
  auto normal_bs0 = std::make_shared<Basis>(normal_bs);
  auto normal_bs1 =
      detail::shift_basis_origin(*normal_bs0, zero_shift_base, R_max, dcell);

  // F = H
  auto F = H;

  // F += (2J - K)
  for (auto RJ = 0; RJ < RJ_size; ++RJ) {
    using ::mpqc::lcao::detail::direct_vector;
    auto vec_RJ = direct_vector(RJ, RJ_max, dcell);
    auto min_bs0 = detail::shift_basis_origin(min_bs, vec_RJ);
    auto min_bs1 = min_bs0;
    // F += 2 J
    DirectTArray &g_J = g_J_vector[RJ];
    if (!g_J.array().is_initialized()) {
      bases = BasisVector{{*normal_bs0, *normal_bs1, *min_bs0, *min_bs1}};
      engine_pool = make_engine_pool(
          libint2::Operator::coulomb,
          utility::make_array_of_refs(bases[0], bases[1], bases[2], bases[3]));
      p_screener =
          detail::make_screener(world, engine_pool, bases, screen, screen_thresh);

      g_J = pao_factory.compute_direct_integrals(world, engine_pool, bases,
                                                 p_screener);
    }
    F("mu, nu") += 2.0 * g_J("mu, nu, lambda, rho") * D("lambda, rho");
    // F -= K
    DirectTArray &g_K = g_K_vector[RJ];
    if (!g_K.array().is_initialized()) {
      bases = BasisVector{{*normal_bs0, *min_bs0, *normal_bs1, *min_bs1}};
      engine_pool = make_engine_pool(
          libint2::Operator::coulomb,
          utility::make_array_of_refs(bases[0], bases[1], bases[2], bases[3]));
      p_screener =
          detail::make_screener(world, engine_pool, bases, screen, screen_thresh);

      g_K = pao_factory.compute_direct_integrals(world, engine_pool, bases,
                                                 p_screener);
    }
    F("mu, nu") -= g_K("mu, lambda, nu, rho") * D("lambda, rho");
  }

  return F;
}

}  // namespace gaussian
}  // namespace lcao
}  // namespace mpqc
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_SOAD_H_
