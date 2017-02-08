//
// Created by Chong Peng on 10/6/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_LCAO_WFN_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_LCAO_WFN_H_

#include "mpqc/chemistry/qc/lcao/expression/orbital_space.h"
#include "mpqc/chemistry/qc/lcao/expression/trange1_engine.h"
#include "mpqc/chemistry/qc/lcao/integrals/lcao_factory.h"
#include "mpqc/chemistry/qc/lcao/integrals/periodic_lcao_factory.h"
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
  using LCAOFactoryType = lcao::LCAOFactory<Tile, Policy>;

  // clang-format off
  /**
   *  \brief The KeyVal constructor
   *
   * The KeyVal object will be queried for all keywords of Wavefunction and
   * LCAOFactory, and the following keywords:
   *
   * | KeyWord | Type | Default| Description |
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
    lcao_factory_ = lcao::detail::construct_lcao_factory<Tile, Policy>(kv);

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

  virtual ~LCAOWavefunction() = default;

  LCAOFactoryType &lcao_factory() { return *lcao_factory_; }
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
    // - free canonical orbs, since these have most information
    // - free populated orbs, since prefer free always
    // - computed populated orbs, since this is less expensive
    // - computed canonical orbs, since this is most expensive
    auto canonical_orbs_provider = std::dynamic_pointer_cast<
        typename CanonicalOrbitalSpace<TArray>::Provider>(ref_wfn);
    auto has_canonical_orbs =
        canonical_orbs_provider &&
        canonical_orbs_provider->can_evaluate();
    auto canonical_orbs_are_free =
        has_canonical_orbs && canonical_orbs_provider->is_available();
    auto populated_orbs_provider = std::dynamic_pointer_cast<
        typename PopulatedOrbitalSpace<TArray>::Provider>(ref_wfn);
    auto has_populated_orbs =
        populated_orbs_provider &&
        populated_orbs_provider->can_evaluate();

    if (canonical_orbs_are_free || (!has_populated_orbs && has_canonical_orbs)) {
      make_closed_shell_canonical_sdref_subspaces(
          lcao_factory_, canonical_orbs_provider, target_ref_precision, ndocc,
          n_frozen_core, occ_block_, unocc_block_);
    }
    // else ask for populated spaces
    else if (has_populated_orbs) {
      make_closed_shell_sdref_subspaces(
          lcao_factory_, populated_orbs_provider, target_ref_precision, ndocc,
          n_frozen_core, occ_block_, unocc_block_);
    } else
      throw ProgrammingError(
          "ref_wfn does not provide canonical or populated orbitals", __FILE__,
          __LINE__);
  }

  const std::shared_ptr<const ::mpqc::utility::TRange1Engine> &trange1_engine()
      const {
    return lcao_factory_->orbital_space().trange1_engine();
  }

  bool is_frozen_core() const { return frozen_core_; }
  int charge() const { return charge_; }
  size_t occ_block() const { return occ_block_; }
  size_t unocc_block() const { return unocc_block_; }

 private:
  std::shared_ptr<LCAOFactoryType> lcao_factory_;
  bool frozen_core_;
  int charge_;
  std::size_t occ_block_;
  std::size_t unocc_block_;
};


#if TA_DEFAULT_POLICY == 0
extern template class LCAOWavefunction<TA::TensorD, TA::DensePolicy>;
#elif TA_DEFAULT_POLICY == 1
extern template class LCAOWavefunction<TA::TensorD, TA::SparsePolicy>;
#endif
}
}

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_LCAO_WFN_H_
