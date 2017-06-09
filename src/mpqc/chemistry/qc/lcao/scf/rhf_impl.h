//
// Created by Chong Peng on 10/5/16.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_SCF_RHF_IMPL_H_
#define SRC_MPQC_CHEMISTRY_QC_SCF_RHF_IMPL_H_

#include "mpqc/chemistry/qc/lcao/scf/rhf.h"

#include <memory>

#include <madness/world/worldmem.h>
#include "mpqc/chemistry/qc/lcao/expression/trange1_engine.h"
#include "mpqc/chemistry/qc/lcao/integrals/integrals.h"
#include "mpqc/chemistry/qc/lcao/scf/cadf_builder.h"
#include "mpqc/chemistry/qc/lcao/scf/diagonalize_for_coeffs.h"
#include "mpqc/chemistry/qc/lcao/scf/eigen_solve_density_builder.h"
#include "mpqc/chemistry/qc/lcao/scf/purification_density_build.h"
#include "mpqc/chemistry/qc/lcao/scf/rij_exact_k_fock_builder.h"
#include "mpqc/chemistry/qc/lcao/scf/soad.h"
#include "mpqc/chemistry/qc/lcao/scf/traditional_df_fock_builder.h"
#include "mpqc/chemistry/qc/lcao/scf/traditional_four_center_fock_builder.h"
#include "mpqc/util/misc/time.h"

namespace mpqc {
namespace lcao {

/**
 *  RHF member functions
 */

template <typename Tile, typename Policy>
RHF<Tile, Policy>::RHF(const KeyVal& kv)
    : AOWavefunction<Tile, Policy>(kv), kv_(kv) {
  auto mol = *this->wfn_world()->atoms();

  // get the molecular charge
  const auto charge = kv.value<int>("charge", 0);
  if (mol.total_atomic_number() <= charge)
    throw InputError(
        "net charge cannot be greater than the total nuclear charge", __FILE__,
        __LINE__, "charge");
  nelectrons_ = mol.total_atomic_number() - charge;
  if (nelectrons_ % 2 != 0)
    throw InputError("RHF requires an even number of electrons", __FILE__,
                     __LINE__, "charge");

  max_iter_ = kv.value<int>("max_iter", 30);

  density_builder_str_ =
      kv.value<std::string>("density_builder", "eigen_solve");
  localize_ = kv.value<bool>("localize", false);
  t_cut_c_ = kv.value<double>("t_cut_c", 0.0);
}

template <typename Tile, typename Policy>
void RHF<Tile, Policy>::init(const KeyVal& kv) {
  auto& ao_factory = this->ao_factory();
  auto& world = this->wfn_world()->world();
  const auto& mol = *this->wfn_world()->atoms();

  // Overlap ints
  S_ = ao_factory.compute(L"<κ|λ>");
  // H core int
  H_ = ao_factory.compute(L"<κ|H|λ>");

  // fock builder
  init_fock_builder();

  // emultipole integral TODO better interface to compute this
  const auto& basis =
      *this->wfn_world()->basis_registry()->retrieve(OrbitalIndex(L"λ"));
  const auto bs_array = utility::make_array(basis, basis);
  auto multi_pool = gaussian::make_engine_pool(
      libint2::Operator::emultipole1, utility::make_array_of_refs(basis));
  auto r_xyz =
      gaussian::xyz_integrals<Tile, Policy>(world, multi_pool, bs_array);

  const auto nocc = nelectrons_ / 2;

  // density builder
  std::size_t n_cluster = mol.nclusters();
  if (density_builder_str_ == "purification") {
    auto density_builder = scf::PurificationDensityBuilder<Tile, Policy>(
        S_, r_xyz, nocc, n_cluster, t_cut_c_, localize_);
    d_builder_ =
        std::make_unique<decltype(density_builder)>(std::move(density_builder));
  } else if (density_builder_str_ == "eigen_solve") {
    std::string decompo_type =
        kv.value<std::string>("decompo_type", "conditioned");
    auto density_builder = scf::ESolveDensityBuilder<Tile, Policy>(
        S_, r_xyz, nocc, n_cluster, t_cut_c_, decompo_type, localize_);
    d_builder_ =
        std::make_unique<decltype(density_builder)>(std::move(density_builder));
  } else {
    throw std::runtime_error("Unknown DensityBuilder name! \n");
  }

  if (!F_.is_initialized()) {
    // soad
    F_ = gaussian::fock_from_soad(world, mol, basis, H_);
  }

  F_diis_ = F_;
  compute_density();
}

template <typename Tile, typename Policy>
void RHF<Tile, Policy>::init_fock_builder() {
  auto& ao_factory = this->ao_factory();
  auto eri4 = ao_factory.compute(L"(μ ν| G|κ λ)");
  auto builder =
      scf::ReferenceFourCenterFockBuilder<Tile, Policy, decltype(eri4)>(eri4,
                                                                        eri4);
  f_builder_ = std::make_unique<decltype(builder)>(std::move(builder));
}

template <typename Tile, typename Policy>
void RHF<Tile, Policy>::obsolete() {
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
  return this->wfn_world()->atoms()->nuclear_repulsion_energy() +
         D_("i,j").dot(F_("i,j") + H_("i,j"), D_.world()).get();
}

template <typename Tile, typename Policy>
void RHF<Tile, Policy>::solve(int64_t max_iters, double thresh) {
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
    throw MaxIterExceeded("RHF SCF did not converge", __FILE__, __LINE__,
                          max_iters);

  // commit the converged state
  this->computed_ = true;
  energy_ = old_energy;
  // store fock matix in registry
  auto& registry = this->ao_factory().registry();
  f_builder_->register_fock(F_, registry);
  // update precision
  computed_precision_ = thresh;
}

template <typename Tile, typename Policy>
void RHF<Tile, Policy>::compute_density() {
  auto dc_pair = d_builder_->operator()(F_diis_);
  D_ = dc_pair.first;
  C_ = dc_pair.second;
}

template <typename Tile, typename Policy>
void RHF<Tile, Policy>::build_F() {
  auto G = f_builder_->operator()(D_, C_);
  F_("i,j") = H_("i,j") + G("i,j");
}

template <typename Tile, typename Policy>
bool RHF<Tile, Policy>::can_evaluate(Energy* energy) {
  // can only evaluate the energy
  return energy->order() == 0;
}

template <typename Tile, typename Policy>
void RHF<Tile, Policy>::do_evaluate(double target_precision) {
  // initialize is computing for the first time
  if (!this->computed()) {
    init(kv_);
  }

  // do SCF if need re-compute due to obsoletion, or need higher precision
  if (!this->computed() || computed_precision_ > target_precision) {
    solve(max_iter_, target_precision);
  }
}

template <typename Tile, typename Policy>
void RHF<Tile, Policy>::evaluate(Energy* result) {
  do_evaluate(result->target_precision(0));
  this->set_value(result, energy_);
}

template <typename Tile, typename Policy>
bool RHF<Tile, Policy>::can_evaluate(CanonicalOrbitalSpace<array_type>*) {
  return density_builder_str_ == "eigen_solve" && !localize_;
}

template <typename Tile, typename Policy>
void RHF<Tile, Policy>::evaluate(CanonicalOrbitalSpace<array_type>* result,
                                 double target_energy_precision,
                                 std::size_t target_blocksize) {
  if (!can_evaluate(result))
    throw ProgrammingError(
        "RHF: canonical orbitals requested, but can't compute them", __FILE__,
        __LINE__);

  do_evaluate(target_energy_precision);

  // Eigen-based density builder provides full set of canonical orbitals, other
  // build them
  auto d_eigen_builder_ =
      dynamic_cast<scf::ESolveDensityBuilder<Tile, Policy>*>(
          &*d_builder_.get());
  assert(d_eigen_builder_ != nullptr && !d_eigen_builder_->localize());

  auto eps = d_eigen_builder_->orbital_energies();
  auto C = d_eigen_builder_->C();

  std::vector<double> eps_vec(eps.rows());
  std::copy(eps.data(), eps.data() + eps.rows(), eps_vec.begin());
  *result = CanonicalOrbitalSpace<array_type>(OrbitalIndex(L"p"),
                                              OrbitalIndex(L"κ"), C, eps_vec);
}

template <typename Tile, typename Policy>
bool RHF<Tile, Policy>::can_evaluate(PopulatedOrbitalSpace<array_type>*) {
  return density_builder_str_ == "eigen_solve";
}

template <typename Tile, typename Policy>
void RHF<Tile, Policy>::evaluate(PopulatedOrbitalSpace<array_type>* result,
                                 double target_energy_precision,
                                 std::size_t target_blocksize) {
  if (!can_evaluate(result))
    throw ProgrammingError(
        "RHF: populatd orbitals requested, but can't compute them", __FILE__,
        __LINE__);

  do_evaluate(target_energy_precision);

  if (!C_.is_initialized()) {
    throw ProgrammingError(
        "RHF: requested PopulatedOrbitalSpace but not produced by the density "
        "builder");
  }

  const auto ndocc = nelectrons_ / 2;
  std::vector<double> occupancies(ndocc, 2.0);
  *result = PopulatedOrbitalSpace<array_type>(
      OrbitalIndex(L"m"), OrbitalIndex(L"κ"), C_, occupancies);
}

/**
 *  RIRHF member functions
 */
template <typename Tile, typename Policy>
RIRHF<Tile, Policy>::RIRHF(const KeyVal& kv) : RHF<Tile, Policy>(kv) {}

template <typename Tile, typename Policy>
void RIRHF<Tile, Policy>::init_fock_builder() {
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
DirectRIRHF<Tile, Policy>::DirectRIRHF(const KeyVal& kv)
    : RHF<Tile, Policy>(kv) {}

template <typename Tile, typename Policy>
void DirectRIRHF<Tile, Policy>::init_fock_builder() {
  auto& ao_factory = this->ao_factory();

  auto inv = ao_factory.compute(L"( Κ | G| Λ )");
  auto eri3 = ao_factory.compute_direct(L"( Κ | G|κ λ)");

  scf::DFFockBuilder<Tile, Policy, decltype(eri3)> builder(inv, eri3);
  this->f_builder_ = std::make_unique<decltype(builder)>(std::move(builder));
}

/**
 * CadfRHF member functions
 */
template <typename Tile, typename Policy>
CadfRHF<Tile, Policy>::CadfRHF(const KeyVal& kv) : RHF<Tile, Policy>(kv) {
  force_shape_threshold_ = kv.value<double>("force_shape_threshold", 0.0);
  auto user_tcutc = kv.exists("tcutc");
  if (user_tcutc) {
    tcutc_ = kv.value<double>("tcutc");
  } else if (force_shape_threshold_ > 0) {
    tcutc_ = 1e-4;
  } else {
    tcutc_ = 0.0;
  }

  secadf_ = kv.value<bool>("secadf", false);
}

template <typename Tile, typename Policy>
void CadfRHF<Tile, Policy>::init_fock_builder() {
  using DirectArray = typename gaussian::AOFactory<Tile, Policy>::DirectTArray;
  using Builder = scf::CADFFockBuilder<Tile, Policy, DirectArray>;
  this->f_builder_ = std::make_unique<Builder>(
      this->ao_factory(), force_shape_threshold_, tcutc_, secadf_);
}

/**
 * DirectRHF member functions
 */
template <typename Tile, typename Policy>
DirectRHF<Tile, Policy>::DirectRHF(const KeyVal& kv) : RHF<Tile, Policy>(kv) {}

template <typename Tile, typename Policy>
void DirectRHF<Tile, Policy>::init_fock_builder() {
  auto& factory = this->ao_factory();
  auto& world = factory.world();
  auto& ao_factory = ::mpqc::lcao::gaussian::to_ao_factory(factory);
  auto screen = ao_factory.screen();
  auto screen_threshold = ao_factory.screen_threshold();
  auto basis =
      this->wfn_world()->basis_registry()->retrieve(OrbitalIndex(L"λ"));
  this->f_builder_ = std::make_unique<scf::FourCenterFockBuilder<Tile, Policy>>(
      world, basis, basis, basis, true, true, screen, screen_threshold);
}

/**
 * DirectRIRHF member functions
 */
template <typename Tile, typename Policy>
RIJEXACTKRHF<Tile, Policy>::RIJEXACTKRHF(const KeyVal& kv)
    : RHF<Tile, Policy>(kv) {}

template <typename Tile, typename Policy>
void RIJEXACTKRHF<Tile, Policy>::init_fock_builder() {
  auto& ao_factory = this->ao_factory();

  auto inv = ao_factory.compute(L"( Κ | G| Λ )");
  auto eri3 = ao_factory.compute_direct(L"( Κ | G|κ λ)[a_bb]");
  auto basis =
      this->wfn_world()->basis_registry()->retrieve(OrbitalIndex(L"λ"));

  using Builder = scf::RIJEXACTKBuilder<Tile, Policy, decltype(eri3)>;
  this->f_builder_ = std::make_unique<Builder>(inv, eri3, basis, basis, basis);
}

}  // namespace lcao
}  // namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_SCF_RHF_IMPL_H_
