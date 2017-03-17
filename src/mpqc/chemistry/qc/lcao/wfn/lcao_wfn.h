//
// Created by Chong Peng on 10/6/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_LCAO_WFN_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_LCAO_WFN_H_

#include "mpqc/chemistry/qc/lcao/expression/orbital_space.h"
#include "mpqc/chemistry/qc/lcao/expression/trange1_engine.h"
#include "mpqc/chemistry/qc/lcao/factory/lcao_factory.h"
#include "mpqc/chemistry/qc/lcao/factory/periodic_lcao_factory.h"
#include "mpqc/chemistry/qc/lcao/wfn/wfn.h"
#include "mpqc/chemistry/qc/properties/property.h"
#include "mpqc/util/keyval/keyval.h"
#include <mpqc/chemistry/qc/lcao/scf/mo_build.h>


namespace mpqc {
namespace lcao {

/// LCAOWavefunction is a Wavefunction with an LCAOFactory

/// This models wave function methods expressed in LCAO basis (e.g. traditional
/// electron correlation methods, like MO-basis CCSD).
/// \todo elaborate LCAOWavefunction documentation
template <typename Tile, typename Policy>
class LCAOWavefunction : public Wavefunction {
 public:
  using ArrayType = TA::DistArray<Tile, Policy>;
  using DirectArrayType = gaussian::DirectArray<Tile, Policy>;
  using LCAOFactoryType = LCAOFactoryBase<Tile,Policy>;
  using AOFactoryType = gaussian::AOFactoryBase<Tile,Policy>;

  // clang-format off
  /**
   *  \brief The KeyVal constructor
   *
   * The KeyVal object will be queried for all keywords of Wavefunction and
   * LCAOFactory, and the following keywords:
   *
   * | Keyword | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | \c "frozen_core" | bool | true | if true, core electrons are not correlated |
   * | \c "charge" | int | 0 | the net charge of the molecule (derived classes may refine the meaning of this keyword) |
   * | \c "obs_block_size" | int | 24 | the target OBS (Orbital Basis Set) space block size |
   * | \c "occ_block_size" | int | \c "$obs_block_size" | the target block size of the occupied space |
   * | \c "unocc_block_size" | int | \c "$obs_block_size" | the target block size of the unoccupied space |
   *
   */
  // clang-format on
  LCAOWavefunction(const KeyVal &kv) : Wavefunction(kv) {
    init_factory(kv);
    frozen_core_ = kv.value<bool>("frozen_core", true);

    charge_ = kv.value<int>("charge", 0);
    if (this->atoms()->total_atomic_number() <= charge_)
      throw InputError(
          "net charge cannot be greater than the sum of atomic numbers",
          __FILE__, __LINE__, "charge");
    const auto nelectrons = this->atoms()->total_atomic_number() - charge_;
    // TODO remove this when have open-shell capability
    if (nelectrons % 2 != 0)
      throw InputError(
          "LCAOWavefunction for now requires an even number of electrons",
          __FILE__, __LINE__, "charge");
    std::size_t mo_block = kv.value<int>("obs_block_size", 24);
    occ_block_ = kv.value<int>("occ_block_size", mo_block);
    unocc_block_ = kv.value<int>("unocc_block_size", mo_block);
  }

  virtual ~LCAOWavefunction() { }

  LCAOFactoryType &lcao_factory() { return *lcao_factory_; }

  const LCAOFactoryType &lcao_factory() const { return *lcao_factory_; }

  AOFactoryType &ao_factory() { return *ao_factory_; }

  const AOFactoryType &ao_factory() const { return *ao_factory_; }

  void obsolete() override {
    lcao_factory_->obsolete();
    Wavefunction::obsolete();
  }

  /// @brief initializes the orbital spaces associated with a (closed-shell, for
  /// now)
  /// single-determinant reference state

  /// @param ref_wfn the reference wave function; expected to provide either
  ///                CanonicalOrbitalSpace or PopulatedOrbitalSpace
  /// @param target_ref_precision the precision with which to compute the
  ///                orbitals (also see OrbitalSpace::Provider::evaluate())
  virtual void init_sdref(std::shared_ptr<lcao::Wavefunction> ref_wfn,
                          double target_ref_precision) {
    using TArray = TA::DistArray<Tile, Policy>;

    std::size_t n_frozen_core = 0;
    if (frozen_core_) {
      auto n_core_electrons = this->atoms()->core_electrons();
      ExEnv::out0() << indent << "Frozen Core: " << n_core_electrons
                    << " electrons" << std::endl;
      n_frozen_core = n_core_electrons / 2;
    }
    const auto nelectrons = this->atoms()->total_atomic_number() - charge_;
    if (nelectrons % 2 != 0)
      throw ProgrammingError(
          "LCAOWavefunction::init_sdref: closed-shell determinant requires an "
          "even number of electrons",
          __FILE__, __LINE__);
    const auto ndocc = nelectrons / 2;

    // dispatch request for either Canonical or Populated orbitals
    // order of preference:
    // - canonical orbs, since these have most information
    // - populated orbs, since don't have to compute anything
    // - compute canonical orbs ourselves, since this costs
    auto canonical_orbs_provider = std::dynamic_pointer_cast<
        typename CanonicalOrbitalSpace<TArray>::Provider>(ref_wfn);
    auto has_canonical_orbs =
        canonical_orbs_provider && canonical_orbs_provider->can_evaluate();
    auto populated_orbs_provider = std::dynamic_pointer_cast<
        typename PopulatedOrbitalSpace<TArray>::Provider>(ref_wfn);
    auto has_populated_orbs =
        populated_orbs_provider && populated_orbs_provider->can_evaluate();

    if (has_canonical_orbs || !has_populated_orbs) {  // use canonical orbs
      auto orbs = std::make_shared<CanonicalOrbitalSpace<TArray>>();
      if (has_canonical_orbs)
        evaluate(*orbs, canonical_orbs_provider, target_ref_precision);
      else  // last resort: compute canonical orbitals
        orbs = make_closed_shell_canonical_orbitals(ao_factory_, ndocc,
                                                    unocc_block_);

      make_closed_shell_canonical_sdref_subspaces(
          lcao_factory_,
          std::const_pointer_cast<const CanonicalOrbitalSpace<TArray>>(orbs),
          ndocc, n_frozen_core, occ_block_, unocc_block_);
    }
    // else ask for populated spaces
    else {
      make_closed_shell_sdref_subspaces(
          ao_factory_, populated_orbs_provider, target_ref_precision, ndocc,
          n_frozen_core, occ_block_, unocc_block_);
    }
  }

  const std::shared_ptr<const ::mpqc::utility::TRange1Engine> &trange1_engine()
      const {
    return lcao_factory_->orbital_registry().trange1_engine();
  }

  bool is_frozen_core() const { return frozen_core_; }
  int charge() const { return charge_; }
  size_t occ_block() const { return occ_block_; }
  size_t unocc_block() const { return unocc_block_; }

 private:
  /**
    *  Default way of initialize factories
    *  use LCAOFactory and AOFactory
    */
  virtual void init_factory(const KeyVal &kv) {
    lcao_factory_ = construct_lcao_factory<Tile, Policy>(kv);
    ao_factory_ = gaussian::construct_ao_factory<Tile, Policy>(kv);
  }

 private:
  std::shared_ptr<LCAOFactoryType> lcao_factory_;
  std::shared_ptr<AOFactoryType> ao_factory_;
  bool frozen_core_;
  int charge_;
  std::size_t occ_block_;
  std::size_t unocc_block_;
};

/// PeriodicLCAOWavefunction is a Wavefunction with a PeriodicLCAOFactory

/// This models wave function methods expressed in LCAO basis (e.g. traditional
/// electron correlation methods, like MO-basis CCSD).
/// \todo elaborate PeriodicLCAOWavefunction documentation
template <typename Tile, typename Policy>
class PeriodicLCAOWavefunction : public Wavefunction {
 public:
  using ArrayType = TA::DistArray<Tile, Policy>;
  using LCAOFactoryType = lcao::PeriodicLCAOFactory<Tile, Policy>;

  // clang-format off
  /**
   *  \brief The KeyVal constructor
   *
   * The KeyVal object will be queried for all keywords of Wavefunction and
   * LCAOFactory, and the following keywords:
   *
   * | Keyword | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | \c "frozen_core" | bool | true | if true, core electrons are not correlated |
   * | \c "obs_block_size" | int | 24 | the target OBS (Orbital Basis Set) space block size |
   * | \c "occ_block_size" | int | \c "$obs_block_size" | the target block size of the occupied space |
   * | \c "unocc_block_size" | int | \c "$obs_block_size" | the target block size of the unoccupied space |
   *
   */
  // clang-format on
  PeriodicLCAOWavefunction(const KeyVal &kv) : Wavefunction(kv) {
    lcao_factory_ =
        lcao::detail::construct_periodic_lcao_factory<Tile, Policy>(kv);

    frozen_core_ = kv.value<bool>("frozen_core", true);
    std::size_t mo_block = kv.value<int>("obs_block_size", 24);
    occ_block_ = kv.value<int>("occ_block_size", mo_block);
    unocc_block_ = kv.value<int>("unocc_block_size", mo_block);
  }

  virtual ~PeriodicLCAOWavefunction() { }

  LCAOFactoryType &lcao_factory() { return *lcao_factory_; }

  void obsolete() override {
    lcao_factory_->obsolete();
    Wavefunction::obsolete();
  }

  const std::shared_ptr<::mpqc::utility::TRange1Engine> trange1_engine() const {
    return trange1_engine_;
  }

  const std::shared_ptr<Eigen::VectorXd> orbital_energy() const {
    return orbital_energy_;
  }

  bool is_frozen_core() const { return frozen_core_; }
  size_t occ_block() const { return occ_block_; }
  size_t unocc_block() const { return unocc_block_; }

 protected:
  std::shared_ptr<Eigen::VectorXd> orbital_energy_;
  std::shared_ptr<::mpqc::utility::TRange1Engine> trange1_engine_;

 private:
  std::shared_ptr<LCAOFactoryType> lcao_factory_;
  bool frozen_core_;
  std::size_t occ_block_;
  std::size_t unocc_block_;
};

#if TA_DEFAULT_POLICY == 0
extern template class LCAOWavefunction<TA::TensorD, TA::DensePolicy>;
#elif TA_DEFAULT_POLICY == 1
extern template class LCAOWavefunction<TA::TensorD, TA::SparsePolicy>;
extern template class PeriodicLCAOWavefunction<TA::TensorZ, TA::SparsePolicy>;
#endif
}
}

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_LCAO_WFN_H_
