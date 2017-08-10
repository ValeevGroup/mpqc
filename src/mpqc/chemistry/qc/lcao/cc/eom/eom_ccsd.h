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

// close-shell eom-ccsd
template <typename Tile, typename Policy>
class EOM_CCSD : public CCSD<Tile, Policy>, public Provides<ExcitationEnergy> {
 public:
  using TArray = TA::DistArray<Tile, Policy>;
  using GuessVector = ::mpqc::cc::T1T2<TArray, TArray>;
  using numeric_type = typename Tile::numeric_type;

 public:
  // clang-format off
  /**
   * KeyVal constructor
   * @param kv
   *
   * | Keyword | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | max_vector | int | 8 | max number of guess vector per root |
   * | vector_threshold | double | 10 * precision of property | threshold for the norm of new guess vector |
   * | eom_pno | string | none | if to simulate pno, avaialble \c default \c state-specific |
   * | eom_pno_canonical | bool | true | if canonicalize PNOs and OSVs |
   * | eom_tpno | double | 0 | PNO truncation threshold for eom |
   * | eom_tosv | double | 0 | OSV truncation threshold for eom |
   *
   */

  // clang-format on
  EOM_CCSD(const KeyVal &kv) : CCSD<Tile, Policy>(kv) {
    max_vector_ = kv.value<int>("max_vector", 8);
    // will be overwrited by precision of property if default set
    vector_threshold_ = kv.value<double>("vector_threshold", 0);
    eom_pno_ = kv.value<std::string>("eom_pno", "");
    if (!eom_pno_.empty() &&
        (eom_pno_ != "default" && eom_pno_ != "state-specific")) {
      throw InputError("Invalid PNO Simulation method in EOM-CCSD! \n",
                       __FILE__, __LINE__, "eom_pno");
    }
    eom_pno_canonical_ = kv.value<bool>("eom_pno_canonical", true);
    eom_tpno_ = kv.value<double>("eom_tpno", 0.0);
    eom_tosv_ = kv.value<double>("eom_tosv", 0.0);
  }

  void obsolete() override {
    CCSD<Tile, Policy>::obsolete();
    TArray g_ijab_ = TArray();

    TArray FAB_ = TArray();
    TArray FIJ_ = TArray();
    TArray FIA_ = TArray();

    TArray WIbAj_ = TArray();
    TArray WIbaJ_ = TArray();

    TArray WAbCd_ = TArray();
    TArray WAbCi_ = TArray();

    TArray WKlIj_ = TArray();
    TArray WKaIj_ = TArray();

    TArray WAkCd_ = TArray();
    TArray WKlIc_ = TArray();
  }

 protected:
  using CCSD<Tile, Policy>::can_evaluate;
  using CCSD<Tile, Policy>::evaluate;

  bool can_evaluate(ExcitationEnergy *ex_energy) override {
    return ex_energy->order() == 0;
  }

  void evaluate(ExcitationEnergy *ex_energy) override;

 private:
  std::size_t max_vector_;   // max number of guess vector
  double vector_threshold_;  // threshold for norm of new guess vector
  std::string eom_pno_;
  bool eom_pno_canonical_;
  double eom_tpno_;
  double eom_tosv_;

  TArray g_ijab_;

  // F intermediates
  TArray FAB_;
  TArray FIJ_;
  TArray FIA_;

  // W intermediates
  TArray WIbAj_;
  TArray WIbaJ_;
  TArray WAbCd_;  // this may not be initialized
  TArray WAbCi_;
  TArray WKlIj_;
  TArray WKaIj_;
  TArray WAkCd_;
  TArray WKlIc_;
  // TArray WKliC_;

  std::vector<GuessVector> C_;  // initial guess vector

  // compute F and W intermediates
  void compute_FWintermediates();

  // compute contractions of HSS, HSD, HDS, and HDD
  //                         with guess vector Ci
  // reference: CPL, 248 (1996), 189
  TArray compute_HSS_HSD_C(const TArray &Cai, const TArray &Cabij);
  TArray compute_HDS_HDD_C(const TArray &Cai, const TArray &Cabij);

  void init() {
    g_ijab_ = this->get_ijab();

    compute_FWintermediates();

    auto remove_integral = [](const Formula &formula) {
      return formula.rank() == 4;
    };

    this->lcao_factory().registry().purge_if(remove_integral);
  }

  EigenVector<numeric_type> eom_ccsd_davidson_solver(std::size_t max_iter,
                                                     double convergence);
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
