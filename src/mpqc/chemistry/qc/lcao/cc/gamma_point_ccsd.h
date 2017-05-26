#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_GAMMA_POINT_CCSD_H
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_GAMMA_POINT_CCSD_H

#include "mpqc/chemistry/qc/lcao/cc/ccsd.h"
#include "mpqc/chemistry/qc/lcao/wfn/ao_wfn.h"

namespace mpqc {
namespace lcao {

template <typename Tile, typename Policy>
class GammaPointCCSD : public CCSD<Tile, Policy> {
 public:
  using LCAOFactoryType = PeriodicLCAOFactory<Tile, Policy>;
  using TArray = typename CCSD<Tile, Policy>::TArray;

  GammaPointCCSD() = default;

  ~GammaPointCCSD() {}

  GammaPointCCSD(const KeyVal &kv) : CCSD<Tile, Policy>(kv), kv_(kv) {
    phf_wfn_ = kv.keyval("ref")
                   .class_ptr<PeriodicAOWavefunction<TA::TensorD, Policy>>();
    lcao_factory_ =
        lcao::detail::construct_periodic_lcao_factory<Tile, Policy>(kv);

    max_iter_ = kv.value<int64_t>("max_iter", 30);
    converge_ = kv.value<double>("converge", 1.0E-7);
    print_detail_ = kv.value<bool>("print_detail", false);
  }

  void obsolete() override {
    gp_ccsd_corr_energy_ = 0.0;
    phf_wfn_->obsolete();
  }

 private:
  std::shared_ptr<PeriodicAOWavefunction<TA::TensorD, Policy>> phf_wfn_;
  std::shared_ptr<LCAOFactoryType> lcao_factory_;

  const KeyVal kv_;
  int64_t max_iter_;
  double converge_;
  bool print_detail_;
  double gp_ccsd_corr_energy_;
  RowMatrixXd C_;
  Eigen::VectorXd eps_;

 private:
  /*!
   * \brief This initializes gamma-point CCSD
   */
  void init_gpccsd() {
    auto nk = phf_wfn_->nk();
    auto k_size = 1 + detail::k_ord_idx(nk(0) - 1, nk(1) - 1, nk(2) - 1, nk);
    auto unitcell = lcao_factory_->pao_factory().unitcell();

    int64_t gamma_point;
    if (k_size % 2 == 1)
      gamma_point = (k_size - 1) / 2;
    else
      throw std::invalid_argument(
          "# of k points must be odd in order to run gamma-point methods");

    C_ = phf_wfn_->co_coeff()[gamma_point].real();
    eps_ = phf_wfn_->co_energy()[gamma_point].real();

    auto trange1_engine = mo_insert_gamma_point(
        *lcao_factory_, C_, unitcell, this->occ_block(), this->unocc_block());
    // TODO try to use LCAOWavefunction::init_sdref
    this->lcao_factory().orbital_registry().set_trange1_engine(*trange1_engine);
    this->set_orbital_energy(eps_);
  }

  bool can_evaluate(Energy *energy) override {
    // can only evaluate the energy
    return energy->order() == 0;
  }

  void evaluate(Energy *result) override {
    if (!this->computed()) {
      /// cast ref_wfn to Energy::Provider
      auto ref_evaluator =
          std::dynamic_pointer_cast<typename Energy::Provider>(phf_wfn_);
      if (ref_evaluator == nullptr) {
        std::ostringstream oss;
        oss << "RefWavefunction in GammaPointCCSD" << phf_wfn_->class_key()
            << " cannot compute Energy" << std::endl;
        throw InputError(oss.str().c_str(), __FILE__, __LINE__);
      }

      ref_evaluator->evaluate(result);

      double ref_energy = this->get_value(result).derivs(0)[0];

      // initialize
      init_gpccsd();

      // set the precision
      if (converge_ == 0.0) {
        // if no user provided converge limit, use the default one from Energy
        this->set_target_precision(result->target_precision(0));
      } else {
        this->set_target_precision(converge_);
      }

      TArray t1, t2;

      gp_ccsd_corr_energy_ = this->compute_ccsd_conventional(t1, t2);

      this->computed_ = true;
      this->set_value(result, ref_energy + gp_ccsd_corr_energy_);
    }
  }

  /// <ab|ij>
  const TArray get_abij() override {
    return lcao_factory_->compute(L"<a b |G|i j>");
  }

  /// <ij|ab>
  const TArray get_ijab() override {
    return lcao_factory_->compute(L"<i j |G|a b>");
  }

  /// <ij|kl>
  const TArray get_ijkl() override {
    return lcao_factory_->compute(L"<i j|G|k l>");
  }

  /// <ab|cd>
  const TArray get_abcd() override {
    return lcao_factory_->compute(L"<a b|G|c d>");
  }

  /// <ia|jb>
  const TArray get_iajb() override {
    return lcao_factory_->compute(L"<i a|G|j b>");
  }

  /// <ia|bc>
  const TArray get_iabc() override {
    return lcao_factory_->compute(L"<i a|G|b c>");
  }

  /// <ai|bc>
  const TArray get_aibc() override {
    return lcao_factory_->compute(L"<a i|G|b c>");
  }

  /// <ij|ak>
  const TArray get_ijak() override {
    return lcao_factory_->compute(L"<i j|G|a k>");
  }

  /// <ij|ka>
  const TArray get_ijka() override {
    return lcao_factory_->compute(L"<i j|G|k a>");
  }

  /// <a|F|i>
  const TArray get_fock_ai() override {
    return lcao_factory_->compute(L"<a|F|i>");
  }
};

#if TA_DEFAULT_POLICY == 0
extern template class GammaPointCCSD<TA::TensorD, TA::DensePolicy>;
#elif TA_DEFAULT_POLICY == 1
extern template class GammaPointCCSD<TA::TensorD, TA::SparsePolicy>;
#endif

}  // namespace lcao
}  // namespace mpqc
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_GAMMA_POINT_CCSD_H
