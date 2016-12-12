//
// Created by Chong Peng on 11/15/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_CCSD_T_F12_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_CCSD_T_F12_H_

#include "mpqc/chemistry/qc/cc/ccsd_t.h"
#include "mpqc/chemistry/qc/f12/ccsd_f12.h"
#include "mpqc/mpqc_config.h"

namespace mpqc {
namespace f12 {

/**
 * \brief CCSD(T)F12 class
 */

template <typename Tile>
class CCSD_T_F12 : public cc::CCSD_T<Tile, TA::SparsePolicy>,
                   public CCSD_F12<Tile> {
 public:
  /**
   * KeyVal constructor
   * @param kv
   *
   * keywords : all keywords from CCSD_T and CCSDF12
   */

  CCSD_T_F12(const KeyVal& kv)
      : cc::CCSD<Tile, TA::SparsePolicy>(kv),
        cc::CCSD_T<Tile, TA::SparsePolicy>(kv),
        CCSD_F12<Tile>(kv) {}

  virtual ~CCSD_T_F12() {}

  double value() override {
    if (this->energy_ == 0.0) {
      auto& world = this->wfn_world()->world();

      // compute CCSD(F12) first
      auto ccsdf12_time0 = mpqc::fenced_now(world);

      double ccsd_f12 = CCSD_F12<Tile>::value();

      auto ccsdf12_time1 = mpqc::fenced_now(world);
      auto ccsdf12_time = mpqc::duration_in_s(ccsdf12_time0, ccsdf12_time1);
      mpqc::utility::print_par(world, "Total CCSD(F12) Time:  ", ccsdf12_time,
                               "\n");

      // compute (T) energy
      this->lcao_factory().ao_factory().registry().purge(world);

      cc::CCSD_T<Tile, TA::SparsePolicy>::compute_ccsd_t();

      auto ccsdtf12_time1 = mpqc::fenced_now(world);
      auto ccsdtf12_time = mpqc::duration_in_s(ccsdf12_time0, ccsdtf12_time1);
      mpqc::utility::print_par(world, "Total CCSD(T)F12 Time:  ", ccsdtf12_time,
                               "\n");

      this->energy_ = ccsd_f12 + this->triples_energy();
    }
    return this->energy_;
  }

  void obsolete() override {
    CCSD_F12<Tile>::obsolete();
    cc::CCSD_T<Tile, TA::SparsePolicy>::obsolete();
  }

  void compute(qc::PropertyBase *pb) override {
    throw std::runtime_error("Not Implemented!!");
  }
};

#if TA_DEFAULT_POLICY == 1
extern template class CCSD_T_F12<TA::TensorD>;
#endif


} //namespace f12
} // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_CCSD_T_F12_H_
