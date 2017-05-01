#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_SOAD_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_SOAD_H_

#include <libint2/chemistry/sto3g_atomic_density.h>
#include <tiledarray.h>

#include "mpqc/chemistry/molecule/unit_cell.h"
#include "mpqc/chemistry/qc/lcao/basis/basis.h"
#include "mpqc/chemistry/qc/lcao/factory/periodic_ao_factory.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_four_center_fock_builder.h"
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
  auto t0 = mpqc::now(world, true);

  using TArray = typename FactoryType::TArray;
  using DirectTArray = typename FactoryType::DirectTArray;
  using Builder = scf::PeriodicFourCenterFockBuilder<Tile, Policy>;

  // soad density
  auto D_eig = soad_density_eig_matrix(unitcell);

  // get minimal basis
  auto min_bs = std::make_shared<const Basis>(
      parallel_make_basis(world, Basis::Factory("sto-3g"), unitcell));

  // transform soad density from Eigen to TA array
  auto trange1 = min_bs->create_trange1();
  auto D =
      array_ops::eigen_to_array<Tile, Policy>(world, D_eig, trange1, trange1);

  // get necessary information for PeriodicFourCenterFockBuilder ctor
  auto dcell = unitcell.dcell();
  auto R_max = pao_factory.R_max();
  auto RJ_max = pao_factory.RJ_max();
  Vector3i RD_max = {0, 0, 0};
  auto R_size = pao_factory.R_size();
  auto RJ_size = pao_factory.RJ_size();
  int64_t RD_size = 1;
  auto screen = pao_factory.screen();
  auto screen_thresh = pao_factory.screen_threshold();

  // get normal basis
  auto normal_bs = pao_factory.basis_registry()->retrieve(OrbitalIndex(L"λ"));

  // F = H
  auto F = H;

  // F += 2J - K
  auto four_center_fock_builder = std::make_unique<Builder>(
      world, normal_bs, min_bs, dcell, R_max, RJ_max, RD_max, R_size, RJ_size,
      RD_size, true, true, screen, screen_thresh);
  auto G = four_center_fock_builder->operator()(
      D, std::numeric_limits<double>::epsilon());
  F("mu, nu") += G("mu, nu");

  auto t1 = mpqc::now(world, true);
  double time = mpqc::duration_in_s(t0, t1);

  ExEnv::out0() << "\nSOAD Time: " << time << " s" << std::endl;

  return F;
}

}  // namespace gaussian
}  // namespace lcao
}  // namespace mpqc
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_SOAD_H_
