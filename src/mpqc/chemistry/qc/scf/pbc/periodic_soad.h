#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_SOAD_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_SOAD_H_

#include <libint2/chemistry/sto3g_atomic_density.h>
#include <tiledarray.h>

#include "mpqc/chemistry/molecule/unit_cell.h"
#include "mpqc/chemistry/qc/basis/basis.h"
#include "mpqc/chemistry/qc/basis/basis_set.h"
#include "mpqc/chemistry/qc/integrals/periodic_ao_factory.h"
#include "mpqc/chemistry/qc/scf/soad.h"

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
 * \param pao_factory PeriodicAOFactory object
 * \param op a functor that takes TA::TensorZ && and returns a valid tile type
 * \return Fock matrix in real space
 */
template <typename Tile, typename Policy, typename FactoryType>
TA::DistArray<Tile,Policy> periodic_fock_soad(
    madness::World &world, UnitCell const &unitcell, TA::DistArray<Tile,Policy> const &H,
    FactoryType &pao_factory,
    std::function<Tile(TA::TensorZ &&)> op = TA::Noop<TA::TensorZ, true>()) {
  if (world.rank() == 0) {
    std::cout << "\nBuilding Fock Matrix from SOAD Guess ...\n";
  }

  auto RJ_size = pao_factory.RJ_size();
  auto RJ_max = pao_factory.RJ_max();
  auto dcell = unitcell.dcell();

  auto F = H;

  // soad density
  auto D_real = soad_density_eig_matrix(unitcell);
  auto D_comp = D_real.cast<std::complex<double>>();

  // get minimal basis
  auto min_bs =
      parallel_construct_basis(world, BasisSet("sto-3g"), unitcell);

  // transform soad density from Eigen to TA
  auto min_bases = BasisVector{{min_bs, min_bs}};
  auto min_trange = detail::create_trange(min_bases);
  auto min_tr0 = min_trange.data()[0];
  auto min_tr1 = min_trange.data()[1];
  auto D = array_ops::eigen_to_array<Tile,Policy>(world, D_comp, min_tr0, min_tr1);

  // get normal basis
  Vector3d zero_shift_base(0.0, 0.0, 0.0);
  auto R_max = pao_factory.R_max();
  auto normal_bs =
      pao_factory.orbital_basis_registry().retrieve(OrbitalIndex(L"λ"));
  auto normal_bs0 = std::make_shared<Basis>(normal_bs);
  auto normal_bs1 =
      detail::shift_basis_origin(*normal_bs0, zero_shift_base, R_max, dcell);

  // F = H + 2J - K
  for (auto RJ = 0; RJ < RJ_size; ++RJ) {
    using ::mpqc::lcao::detail::direct_vector;
    auto vec_RJ = direct_vector(RJ, RJ_max, dcell);
    auto min_bs0 = detail::shift_basis_origin(min_bs, vec_RJ);
    auto min_bs1 = min_bs0;
    // F += 2 J
    auto bases =
        BasisVector{{*normal_bs0, *normal_bs1, *min_bs0, *min_bs1}};
    auto eri_e = make_engine_pool(
        libint2::Operator::coulomb,
        utility::make_array_of_refs(bases[0], bases[1], bases[2], bases[3]));
    auto J = pao_factory.compute_integrals(world, eri_e, bases);
    F("mu, nu") += 2.0 * J("mu, nu, lambda, rho") * D("lambda, rho");
    // F -= K
    auto bases_K =
        BasisVector{{*normal_bs0, *min_bs0, *normal_bs1, *min_bs1}};
    auto eri_e_K = make_engine_pool(
        libint2::Operator::coulomb,
        utility::make_array_of_refs(bases[0], bases[1], bases[2], bases[3]));
    auto K = pao_factory.compute_integrals(world, eri_e_K, bases_K);
    F("mu, nu") -= K("mu, lambda, nu, rho") * D("lambda, rho");
  }

  return F;
}

}  // namespace gaussian
}  // namespace lcao
}  // namespace mpqc
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_SOAD_H_
