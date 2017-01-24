//
// Created by Chong Peng on 11/15/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_CCSD_T_F12_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_CCSD_T_F12_H_

#include "mpqc/chemistry/qc/cc/ccsd_t.h"
#include "mpqc/chemistry/qc/f12/ccsd_f12.h"
#include "mpqc/mpqc_config.h"

namespace mpqc {
namespace lcao {

/**
 * \brief CCSD(T)F12 class
 *
 * KeyVal type: CCSD(T)F12
 */

template <typename Tile>
class CCSD_T_F12 : public CCSD_T<Tile, TA::SparsePolicy>,
                   public CCSD_F12<Tile> {
 public:

  // clang-format off
  /**
   * KeyVal constructor
   * @param kv
   *
   * KeyVal keywords : inherit all keywords from CCSD_T and CCSDF12
   */
  // clang-format on

  CCSD_T_F12(const KeyVal& kv)
      : CCSD<Tile, TA::SparsePolicy>(kv),
        CCSD_T<Tile, TA::SparsePolicy>(kv),
        CCSD_F12<Tile>(kv) {}

  virtual ~CCSD_T_F12() {}


  void obsolete() override {
    CCSD_F12<Tile>::obsolete();
    CCSD_T<Tile, TA::SparsePolicy>::obsolete();
  }

 protected:

  void evaluate(Energy* result) override {
    if (!this->computed()) {
      auto& world = this->wfn_world()->world();

      // compute CCSD(F12) first
      CCSD_F12<Tile>::evaluate(result);
      double ccsd_f12_energy = this->get_value(result).derivs(0)[0];

      auto t_time0 = mpqc::fenced_now(world);
      // compute (T) energy
      this->lcao_factory().ao_factory().registry().purge(world);
      CCSD_T<Tile, TA::SparsePolicy>::compute_ccsd_t();

      auto t_time1 = mpqc::fenced_now(world);
      auto t_time = mpqc::duration_in_s(t_time0, t_time1);
      mpqc::utility::print_par(world, "(T) Time in CCSD(T)F12:  ", t_time, "\n");

      this->computed_ = true;
      this->set_value(result, ccsd_f12_energy + this->triples_energy());
    }
  }

};

#if TA_DEFAULT_POLICY == 1
extern template class CCSD_T_F12<TA::TensorD>;
#endif


}  //namespace lcao
}  // namespace  mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_CCSD_T_F12_H_
