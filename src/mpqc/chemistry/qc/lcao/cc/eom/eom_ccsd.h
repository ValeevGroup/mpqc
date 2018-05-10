//
//
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_CC_EOM_EOM_CCSD_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_CC_EOM_EOM_CCSD_H_

#include "mpqc/chemistry/qc/lcao/cc/ccsd.h"
#include "mpqc/chemistry/qc/lcao/cc/ccsd_hbar.h"
#include "mpqc/chemistry/qc/lcao/cc/eom/eom_preconditioner.h"
#include "mpqc/chemistry/qc/lcao/ci/cis.h"
#include "mpqc/chemistry/qc/properties/excitation_energy.h"

namespace mpqc {
namespace lcao {

/// Implements closed-shell EOM-EA-CCSD
template <typename Tile, typename Policy>
class EOM_CCSD : public CCSD<Tile, Policy>, public Provides<ExcitationEnergy> {
 public:
  using TArray = TA::DistArray<Tile, Policy>;
  using GuessVector = ::mpqc::cc::TPack<TArray>;
  using numeric_type = typename Tile::numeric_type;

 public:
  // clang-format off
  /**
   * The KeyVal constructor.
   *
   * @param kv KeyVal object; it will be queried for all keywords of ::mpqc::lcao::CCSD , as well as the following additional keywords:
   *
   * | Keyword | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | davidson_solver | string | multi-state | choose the davidson solver to use, multi-state or single-state  |
   * | max_vector | int | 8 | max number of guess vector per root |
   * | vector_threshold | real | 10 * precision of property | threshold for the norm of new guess vector |
   * | eom_pno | string | none | if to simulate pno, available \c default, which uses first excited state to generate PNOs \c state-average, use average of states to generate PNOs |
   * | eom_pno_canonical | bool | true | if canonicalize PNOs and OSVs |
   * | eom_tpno | real | 0 | PNO truncation threshold for eom |
   * | eom_tosv | real | 0 | OSV truncation threshold for eom |
   *
   */

  // clang-format on
  EOM_CCSD(const KeyVal &kv) : CCSD<Tile, Policy>(kv) {
    max_vector_ = kv.value<int>("max_vector", 8);
    // will be overwrited by precision of property if default set
    vector_threshold_ = kv.value<double>("vector_threshold", 0);

    // need to modify the method keyword, CIS only has standard and df
    KeyVal &kv_nonconst = const_cast<KeyVal &>(kv);
    std::string original_method = kv.value<std::string>("method", "");
    std::string cis_method = (this->df_ ? "df" : "standard");
    kv_nonconst.assign("method", cis_method);

    // construct CIS Wavefunction
    cis_guess_wfn_ = std::make_shared<CIS<Tile, Policy>>(kv);

    // change method keyword back to original value
    if (!original_method.empty()) {
      kv_nonconst.assign("method", original_method);
    }

    davidson_solver_ = kv.value<std::string>("davidson_solver", "multi-state");

    if (davidson_solver_ != "multi-state" &&
        davidson_solver_ != "single-state") {
      throw InputError("Invalid Davidson Solver in EOM-CCSD! \n", __FILE__,
                       __LINE__, "davidson_solver");
    }

    eom_pno_ = kv.value<std::string>("eom_pno", "");
    if (!eom_pno_.empty() &&
        (eom_pno_ != "default" && eom_pno_ != "state-average")) {
      throw InputError("Invalid PNO Simulation method in EOM-CCSD! \n",
                       __FILE__, __LINE__, "eom_pno");
    }
    eom_pno_canonical_ = kv.value<bool>("eom_pno_canonical", true);
    eom_tpno_ = kv.value<double>("eom_tpno", 0.0);
    eom_tosv_ = kv.value<double>("eom_tosv", 0.0);
  }

  void obsolete() override {
    CCSD<Tile, Policy>::obsolete();
    cis_guess_wfn_->obsolete();
  }

 protected:
  using CCSD<Tile, Policy>::can_evaluate;
  using CCSD<Tile, Policy>::evaluate;

  bool can_evaluate(ExcitationEnergy *ex_energy) override {
    return ex_energy->order() == 0;
  }

  void evaluate(ExcitationEnergy *ex_energy) override;

 private:
  EigenVector<numeric_type> eom_ccsd_davidson_solver(
      std::size_t n_roots, const std::vector<TArray> &cis_vector,
      const std::vector<numeric_type> &cis_eigs, std::size_t max_iter,
      double convergence);
  // compute F and W intermediates
  cc::Intermediates<TArray> compute_FWintermediates();

  // compute contractions of HSS, HSD, HDS, and HDD
  //                         with guess vector Ci
  // reference: CPL, 248 (1996), 189
  TArray compute_HSS_HSD_C(const TArray &Cai, const TArray &Cabij,
                           const cc::Intermediates<TArray> &imds);
  TArray compute_HDS_HDD_C(const TArray &Cai, const TArray &Cabij,
                           const cc::Intermediates<TArray> &imds);

  cc::Intermediates<TArray> init() {
    auto imds = compute_FWintermediates();

    auto remove_integral = [](const Formula &formula) {
      return formula.rank() == 4;
    };

    this->lcao_factory().registry().purge_if(remove_integral);

    return imds;
  }

 private:
  std::shared_ptr<CIS<Tile, Policy>> cis_guess_wfn_;
  std::size_t max_vector_;   // max number of guess vector
  double vector_threshold_;  // threshold for norm of new guess vector
  std::string davidson_solver_;
  std::string eom_pno_;
  bool eom_pno_canonical_;
  double eom_tpno_;
  double eom_tosv_;
};

#if TA_DEFAULT_POLICY == 0
extern template class EOM_CCSD<TA::TensorD, TA::DensePolicy>;
#elif TA_DEFAULT_POLICY == 1
extern template class EOM_CCSD<TA::TensorD, TA::SparsePolicy>;
#endif

}  // namespace lcao
}  // namespace mpqc

#include "eom_ccsd_impl.h"

#endif
