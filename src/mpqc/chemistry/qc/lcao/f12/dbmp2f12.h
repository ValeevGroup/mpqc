//
// Created by Chong Peng on 7/11/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_DBMP2F12_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_DBMP2F12_H_

#include "mpqc/chemistry/qc/lcao/f12/db_f12_intermediates.h"
#include "mpqc/chemistry/qc/lcao/f12/mp2f12.h"
#include "mpqc/chemistry/qc/lcao/mbpt/dbmp2.h"
#include "mpqc/mpqc_config.h"

namespace mpqc {
namespace lcao {


/**
 * \brief Dual Basis MP2F12 method for closed shell with RI
 */

template <typename Tile>
class RIDBRMP2F12 : public RIRMP2F12<Tile> {
 public:
  using TArray = TA::DistArray<Tile,TA::SparsePolicy>;
  /**
   * KeyVal constructor
   * @param kv
   *
   * keywords:  takes all keywords from RIRMP2F12 class
   *
   * | Keyword | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | mp2 | string | none | if to recompute mp2 energy(redo) or update  |
   */
  RIDBRMP2F12(const KeyVal& kv);
  virtual ~RIDBRMP2F12() = default;

  /// override the evaluate function in RMP2F12
  void evaluate(Energy* result) override;

 private:

  /// override the initialize function in RMP2F12
  void init() override;

  TArray compute_B() override;
  TArray compute_V() override;
  TArray compute_X() override;
  double compute_cabs_singles() override;
  double compute_new_mp2();

 private:
  const KeyVal kv_;
  std::string mp2_method_;
};

#if TA_DEFAULT_POLICY == 1
extern template class RIDBRMP2F12<TA::TensorD>;
#endif

}  // namespace lcao
}  // namespace mpqc

#include "dbmp2f12_impl.h"

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_DBMP2F12_H_
