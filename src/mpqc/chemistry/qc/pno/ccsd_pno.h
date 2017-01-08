/*
 * ccsd_pno.h
 *
 *  Created on: Jan 4, 2017
 *      Author: jinmei
 */

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_PNO_CCSD_PNO_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_PNO_CCSD_PNO_H_

#include <tiledarray.h>

#include "mpqc/chemistry/qc/scf/mo_build.h"
#include "mpqc/chemistry/qc/mbpt/denom.h"
#include "mpqc/chemistry/qc/wfn/lcao_wfn.h"

namespace mpqc {
namespace lcao {

  template <typename Tile, typename Policy>
  class CCSD_PNO : public LCAOWavefunction<Tile, Policy> {

  private:
   //const KeyVal kv_;
   std::shared_ptr<Wavefunction> ref_wfn_;
   bool df_;

   void init();

   TA::DistArray<Tile, Policy> compute_mp2_t2();

  public:
   CCSD_PNO() = default;
   CCSD_PNO(const KeyVal &kv);

   ~CCSD_PNO() = default;
   void compute(PropertyBase* pb) override;
   // compute function
   double value() override;

  };
}  // namespace lcao
}  // namespace mpqc

#include "ccsd_pno_impl.h"

#endif // MPQC4_SRC_MPQC_CHEMISTRY_QC_PNO_CCSD_PNO_H_
