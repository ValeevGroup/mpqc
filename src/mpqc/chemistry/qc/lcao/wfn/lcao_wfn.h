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

    // initialize the orbital spaces
    std::shared_ptr<CanonicalOrbitalSpace<TArray>> obs_space =
        std::make_shared<CanonicalOrbitalSpace<TArray>>();
    evaluate(*obs_space, ref_wfn, target_ref_precision);

    auto &lcao_factory = this->lcao_factory();
    auto &orbital_registry = lcao_factory.orbital_space();
    auto &world = lcao_factory.world();

    auto mo_time0 = mpqc::fenced_now(world);

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

    // divide the LCAO space into subspaces using Eigen .. boo
    RowMatrixXd C_obs = array_ops::array_to_eigen(obs_space->coefs());
    const auto nobs = C_obs.cols();
    const auto nao = C_obs.rows();
    RowMatrixXd C_occ = C_obs.block(0, 0, C_obs.rows(), ndocc);
    RowMatrixXd C_corr_occ =
        C_occ.block(0, n_frozen_core, nao, ndocc - n_frozen_core);
    RowMatrixXd C_unocc = C_obs.rightCols(nobs - ndocc);

    ExEnv::out0() << indent << "OccBlockSize: " << occ_block_ << std::endl;
    ExEnv::out0() << indent << "UnoccBlockSize: " << unocc_block_ << std::endl;

    using TRange1Engine = ::mpqc::utility::TRange1Engine;
    auto tre = std::make_shared<TRange1Engine>(ndocc, nobs, occ_block_,
                                               unocc_block_, n_frozen_core);

    // get all the trange1s
    auto tr_ao = obs_space->coefs().trange().data()[0];
    auto tr_corr_occ = tre->get_active_occ_tr1();
    auto tr_occ = tre->compute_range(ndocc, occ_block_);
    auto tr_vir = tre->get_vir_tr1();
    auto tr_all = tre->get_all_tr1();

    mpqc::detail::parallel_print_range_info(world, tr_occ, "Occ");
    mpqc::detail::parallel_print_range_info(world, tr_corr_occ, "CorrOcc");
    mpqc::detail::parallel_print_range_info(world, tr_vir, "Vir");
    mpqc::detail::parallel_print_range_info(world, tr_all, "Obs");

    // convert eigen arrays to TA
    auto C_occ_ta =
        array_ops::eigen_to_array<Tile, Policy>(world, C_occ, tr_ao, tr_occ);
    auto C_corr_occ_ta = array_ops::eigen_to_array<Tile, Policy>(
        world, C_corr_occ, tr_ao, tr_corr_occ);
    auto C_unocc_ta =
        array_ops::eigen_to_array<Tile, Policy>(world, C_unocc, tr_ao, tr_vir);

    // insert orbital spaces into registry
    using OrbitalSpaceTArray = OrbitalSpace<TA::DistArray<Tile, Policy>>;
    auto occ_space =
        OrbitalSpaceTArray(OrbitalIndex(L"m"), OrbitalIndex(L"κ"), C_occ_ta);
    orbital_registry.add(occ_space);

    auto corr_occ_space = OrbitalSpaceTArray(OrbitalIndex(L"i"),
                                             OrbitalIndex(L"κ"), C_corr_occ_ta);
    orbital_registry.add(corr_occ_space);

    auto vir_space =
        OrbitalSpaceTArray(OrbitalIndex(L"a"), OrbitalIndex(L"κ"), C_unocc_ta);
    orbital_registry.add(vir_space);

    orbital_registry.add(*obs_space);

    auto mo_time1 = mpqc::fenced_now(world);
    auto mo_time = mpqc::duration_in_s(mo_time0, mo_time1);
    utility::print_par(world, "closed-shell OBS MO Build Time: ", mo_time,
                       " S\n");
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
   * | KeyWord | Type | Default| Description |
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

  virtual ~PeriodicLCAOWavefunction() = default;

  LCAOFactoryType &lcao_factory() { return *lcao_factory_; }
  void obsolete() override {
    lcao_factory_->registry().purge();
    lcao_factory_->orbital_space().clear();
    lcao_factory_->ao_factory().registry().purge();
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
