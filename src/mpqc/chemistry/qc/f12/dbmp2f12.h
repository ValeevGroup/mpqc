//
// Created by Chong Peng on 7/11/16.
//

#ifndef MPQC_DBMP2F12_H
#define MPQC_DBMP2F12_H

#include <mpqc/chemistry/qc/f12/db_f12_intermediates.h>
#include <mpqc/chemistry/qc/f12/mp2f12.h>
#include <mpqc/chemistry/qc/mbpt/dbmp2.h>

namespace mpqc {
namespace f12 {

class RIDBRMP2F12 : public RIRMP2F12{

public:

  /**
   * KeyVal constructor
   * @param kv
   *
   * keywords:  takes all keywords from RIRMP2F12 class
   *
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | redo_mp2 | bool | false | if to recompute mp2 energy in VBS |
   */
  RIDBRMP2F12(const KeyVal& kv);
  ~RIDBRMP2F12() = default;

  double value() override;

private:

  TArray compute_B() override ;
  TArray compute_V() override ;
  TArray compute_X() override ;
  double compute_cabs_singles() override;

private:

  const KeyVal kv_;
  bool redo_mp2_;

};
}  // namespace f12
}  // namespace mpqc

#endif  // MPQC_DBMP2F12_H
