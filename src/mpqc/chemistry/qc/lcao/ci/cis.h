//
// Created by Chong Peng on 3/1/17.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_LCAO_CI_CIS_H_
#define SRC_MPQC_CHEMISTRY_QC_LCAO_CI_CIS_H_

#include "mpqc/chemistry/qc/lcao/wfn/lcao_wfn.h"
#include "mpqc/chemistry/qc/properties/excitation_energy.h"
#include "mpqc/mpqc_config.h"

namespace mpqc {
namespace lcao {

/**
 * CIS for closed shell system
 *
 */
template <typename Tile, typename Policy>
class CIS : public LCAOWavefunction<Tile, Policy>,
            public Provides<ExcitationEnergy> {
 public:
  /**
  * KeyVal constructor
  * @param kv
  * | Keyword | Type | Default| Description |
  * |---------|------|--------|-------------|
  * | ref | Wavefunction | none | reference Wavefunction, RHF for example |
  */
  explicit CIS(const KeyVal& kv) : LCAOWavefunction<Tile, Policy>(kv) {
    if (kv.exists("ref")) {
      ref_wfn_ = kv.class_ptr<Wavefunction>("ref");
    } else {
      throw InputError("Default Ref Wfn in CIS is not support! \n", __FILE__,
                       __LINE__, "ref");
    }
  }

  ~CIS() = default;

 private:
  bool can_evaluate(ExcitationEnergy* ex_energy) override {
    return ex_energy->order() == 0;
  }

  void evaluate(ExcitationEnergy* ex_energy) override {
    if (!this->computed()) {
      auto n_roots = ex_energy->n_roots();

      std::vector<double> result;
      for (auto i = 0; i < n_roots; i++) {
        result.push_back(i);
      }

      this->computed_ = true;
      this->set_value(ex_energy, result);
    }
  }

 private:
  std::shared_ptr<Wavefunction> ref_wfn_;
};

}  // namespace lcao
}  // namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_LCAO_CI_CIS_H_
