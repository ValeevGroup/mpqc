//
// Created by Chong Peng on 10/5/16.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_SCF_RHF_IMPL_H_
#define SRC_MPQC_CHEMISTRY_QC_SCF_RHF_IMPL_H_

#include "mpqc/chemistry/qc/lcao/integrals/integrals.h"
#include "mpqc/chemistry/qc/lcao/scf/diagonalize_for_coeffs.h"
#include "mpqc/chemistry/qc/lcao/scf/traditional_df_fock_builder.h"
#include "mpqc/chemistry/qc/lcao/scf/traditional_four_center_fock_builder.h"
#include "mpqc/util/external/c++/memory"
#include "mpqc/util/misc/time.h"
#include "purification_density_build.h"
#include "eigen_solve_density_builder.h"
#include "soad.h"
#include <madness/world/worldmem.h>

namespace mpqc {
namespace lcao {

/**
 *  RHF member functions
 */

template <typename Tile, typename Policy>
RHF<Tile,Policy>::RHF(const KeyVal& kv) : AOWavefunction<Tile,Policy>(kv), kv_(kv) {}

template <typename Tile, typename Policy>
void RHF<Tile,Policy>::init(const KeyVal& kv) {
  auto& ao_factory = this->ao_factory();
  auto& world = ao_factory.world();
  auto& mol = ao_factory.molecule();

  // check if even number of electron first
  std::size_t occ = mol.occupation();
  if (occ % 2 != 0) {
    throw std::invalid_argument(
        "RHF doesn't support ODD number of electrons! \n");
  }
  occ = occ / 2;

  max_iter_ = kv.value<int>("max_iter", 30);

  // Overlap ints
  S_ = ao_factory.compute(L"<κ|λ>");
  // H core int
  H_ = ao_factory.compute(L"<κ|H|λ>");

  // fock builder
  init_fock_builder();

  // emultipole integral TODO better interface to compute this
  auto basis = *ao_factory.orbital_basis_registry().retrieve(OrbitalIndex(L"λ"));
  const auto bs_array = utility::make_array(basis, basis);
  auto multi_pool = gaussian::make_engine_pool(
      libint2::Operator::emultipole1, utility::make_array_of_refs(basis));
  auto r_xyz = gaussian::xyz_integrals<Tile, Policy>(world, multi_pool, bs_array);

  // density builder
  std::string density_builder =
      kv.value<std::string>("density_builder", "eigen_solve");
  bool localize = kv.value<bool>("localize", false);
  double t_cut_c = kv.value<double>("t_cut_c", 0.0);
  std::size_t n_cluster = mol.nclusters();
  if (density_builder == "purification") {
    auto density_builder = scf::PurificationDensityBuilder<Tile,Policy>(
        S_, r_xyz, occ, n_cluster, t_cut_c, localize);
    d_builder_ =
        std::make_unique<decltype(density_builder)>(std::move(density_builder));
  } else if (density_builder == "eigen_solve") {
    std::string decompo_type =
        kv.value<std::string>("decompo_type", "conditioned");
    auto density_builder = scf::ESolveDensityBuilder<Tile,Policy>(
        S_, r_xyz, occ, n_cluster, t_cut_c, decompo_type, localize);
    d_builder_ =
        std::make_unique<decltype(density_builder)>(std::move(density_builder));
  } else {
    throw std::runtime_error("Unknown DensityBuilder name! \n");
  }

  if(!F_.is_initialized()){
    // soad
    auto eri_e = gaussian::make_engine_pool(libint2::Operator::coulomb,
                                             utility::make_array_of_refs(basis));
    F_ = gaussian::fock_from_soad(world, mol, basis, eri_e, H_);
  }

  F_diis_ = F_;
  compute_density();
}

template <typename Tile, typename Policy>
void RHF<Tile,Policy>::init_fock_builder() {
  auto& ao_factory = this->ao_factory();
  auto eri4 = ao_factory.compute(L"(μ ν| G|κ λ)");
  auto builder = scf::FourCenterBuilder<Tile, Policy, decltype(eri4)>(std::move(eri4));
  f_builder_ = std::make_unique<decltype(builder)>(std::move(builder));
}

template <typename Tile, typename Policy>
void RHF<Tile,Policy>::obsolete() {
  ::mpqc::Wavefunction::obsolete();

  H_ = array_type();
  S_ = array_type();
  F_ = array_type();
  F_diis_ = array_type();
  D_ = array_type();
  C_ = array_type();

  rhf_times_ = std::vector<double>();
  d_times_ = std::vector<double>();
  build_times_ = std::vector<double>();

  AOWavefunction<Tile, Policy>::obsolete();
}

template <typename Tile, typename Policy>
double RHF<Tile, Policy>::compute_energy() const {
  return this->ao_factory().molecule().nuclear_repulsion() +
         D_("i,j").dot(F_("i,j") + H_("i,j"), D_.world()).get();
}

template <typename Tile, typename Policy>
void RHF<Tile,Policy>::solve(int64_t max_iters, double thresh) {

  // for now use same precision for energy and orbital gradient
  const auto target_energy_precision = thresh;
  const auto target_orbgrad_precision = thresh;

  auto& world = F_.world();

  if (world.rank() == 0) {
    std::cout << "Starting SCF:\n"
              << "\tThreshold: " << thresh << "\n"
              << "\tMaximum number of iterations: " << max_iters << "\n";
  }

  madness::print_meminfo(world.rank(), "RHF:begin");

  TiledArray::DIIS<array_type> diis;

  auto iter = 0;
  auto error = std::numeric_limits<double>::max();
  auto rms_error = std::numeric_limits<double>::max();
  auto old_energy = compute_energy();
  const double volume = F_.trange().elements_range().volume();

  while (iter < max_iters &&
         (target_energy_precision < (error / old_energy) ||
          target_orbgrad_precision < (rms_error / volume))) {
    auto s0 = mpqc::fenced_now(world);

    madness::print_meminfo(world.rank(), "RHF:before_fock");
    build_F();
    auto b1 = mpqc::fenced_now(world);
    madness::print_meminfo(world.rank(), "RHF:after_fock");

    auto current_energy = compute_energy();
    error = std::abs(old_energy - current_energy);
    old_energy = current_energy;

    array_type Grad;
    Grad("i,j") =
        F_("i,k") * D_("k,l") * S_("l,j") - S_("i,k") * D_("k,l") * F_("l,j");
    madness::print_meminfo(world.rank(), "RHF:orbgrad");

    rms_error = Grad("i,j").norm().get();

    F_diis_ = F_;
    diis.extrapolate(F_diis_, Grad);
    madness::print_meminfo(world.rank(), "RHF:diis");

    auto d0 = mpqc::fenced_now(world);
    compute_density();
    auto s1 = mpqc::fenced_now(world);
    madness::print_meminfo(world.rank(), "RHF:density");

    rhf_times_.push_back(mpqc::duration_in_s(s0, s1));
    d_times_.push_back(mpqc::duration_in_s(d0, s1));
    build_times_.push_back(mpqc::duration_in_s(s0, b1));

    if (world.rank() == 0) {
      std::cout << "iteration: " << iter << "\n"
                << "\tEnergy: " << old_energy << "\n"
                << "\tabs(Energy Change)/energy: "
                << (error / std::abs(old_energy)) << "\n"
                << "\t(Gradient Norm)/n^2: " << (rms_error / volume) << "\n"
                << "\tScf Time: " << rhf_times_.back() << "\n"
                << "\t\tDensity Time: " << d_times_.back() << "\n"
                << "\t\tFock Build Time: " << build_times_.back() << "\n";
    }
    f_builder_->print_iter("\t\t");
    d_builder_->print_iter("\t\t");
    ++iter;
  }

  if (iter == max_iters)
    throw MaxIterExceeded("RHF SCF did not converge", __FILE__, __LINE__, max_iters);

  energy_ = old_energy;
  // store fock matix in registry
  auto& registry = this->ao_factory().registry();
  f_builder_->register_fock(F_, registry);
}

template <typename Tile, typename Policy>
void RHF<Tile,Policy>::compute_density() {
  auto dc_pair = d_builder_->operator()(F_diis_);
  D_ = dc_pair.first;
  C_ = dc_pair.second;
}

template <typename Tile, typename Policy>
void RHF<Tile,Policy>::build_F() {
  auto G = f_builder_->operator()(D_, C_);
  F_("i,j") = H_("i,j") + G("i,j");
}

template <typename Tile, typename Policy>
bool RHF<Tile,Policy>::can_evaluate(Energy* energy) {
  // can only evaluate the energy
  return energy->order() == 0;
}

template <typename Tile, typename Policy>
void RHF<Tile,Policy>::evaluate(Energy* result) {
  if(!this->computed()){
    init(kv_);
    solve(max_iter_, result->target_precision(0));
    this->computed_ = true;
    set_value(result, energy_);
  }
}

/**
 *  RIRHF member functions
 */
template <typename Tile, typename Policy>
RIRHF<Tile,Policy>::RIRHF(const KeyVal& kv) : RHF<Tile,Policy>(kv) {}

template <typename Tile, typename Policy>
void RIRHF<Tile,Policy>::init_fock_builder() {
  auto& ao_factory = this->ao_factory();
  auto inv = ao_factory.compute(L"( Κ | G| Λ )");
  auto eri3 = ao_factory.compute(L"( Κ | G|κ λ)");
  scf::DFFockBuilder<Tile, Policy, decltype(eri3)> builder(inv, eri3);
  this->f_builder_ = std::make_unique<decltype(builder)>(std::move(builder));
}

/**
 * DirectRIRHF member functions
 */
template <typename Tile, typename Policy>
DirectRIRHF<Tile,Policy>::DirectRIRHF(const KeyVal& kv) : RHF<Tile,Policy>(kv) {}

template <typename Tile, typename Policy>
void DirectRIRHF<Tile,Policy>::init_fock_builder() {
  auto& direct_ao_factory = this->direct_ao_factory();
  auto& ao_factory = this->ao_factory();

  auto inv = ao_factory.compute(L"( Κ | G| Λ )");
  auto eri3 = direct_ao_factory.compute(L"( Κ | G|κ λ)");

  scf::DFFockBuilder<Tile, Policy, decltype(eri3)> builder(inv, eri3);
  this->f_builder_ = std::make_unique<decltype(builder)>(std::move(builder));
}

/**
 * DirectRHF member functions
 */
template <typename Tile, typename Policy>
DirectRHF<Tile,Policy>::DirectRHF(const KeyVal& kv) : RHF<Tile,Policy>(kv) {}

template <typename Tile, typename Policy>
void DirectRHF<Tile,Policy>::init_fock_builder() {
  auto& direct_ao_factory = this->direct_ao_factory();
  auto eri4 = direct_ao_factory.compute(L"(μ ν| G|κ λ)");
  auto builder = scf::FourCenterBuilder<Tile, Policy, decltype(eri4)>(std::move(eri4));
  this->f_builder_ = std::make_unique<decltype(builder)>(std::move(builder));
}

}  // namespace lcao
}  // namespace mpqc

#endif //SRC_MPQC_CHEMISTRY_QC_SCF_RHF_IMPL_H_
