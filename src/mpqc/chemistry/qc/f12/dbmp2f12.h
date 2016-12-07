//
// Created by Chong Peng on 7/11/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_DBMP2F12_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_DBMP2F12_H_

#include "mpqc/chemistry/qc/f12/db_f12_intermediates.h"
#include "mpqc/chemistry/qc/f12/mp2f12.h"
#include "mpqc/chemistry/qc/mbpt/dbmp2.h"

namespace mpqc {
namespace f12 {


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
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | mp2 | string | none | if to recompute mp2 energy(redo) or update  |
   */
  RIDBRMP2F12(const KeyVal& kv);
  virtual ~RIDBRMP2F12() = default;

  double value() override;

 private:
  TArray compute_B() override;
  TArray compute_V() override;
  TArray compute_X() override;
  double compute_cabs_singles() override;
  double compute_new_mp2();

 private:
  const KeyVal kv_;
  std::string mp2_method_;
};

extern template class RIDBRMP2F12<TA::TensorD>;

}  // namespace f12
}  // namespace mpqc

#include "dbmp2f12_impl.h"

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_DBMP2F12_H_
