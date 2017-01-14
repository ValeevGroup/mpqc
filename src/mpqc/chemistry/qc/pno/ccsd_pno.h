/*
 * ccsd_pno.h
 *
 *  Created on: Jan 4, 2017
 *      Author: jinmei
 */

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_PNO_CCSD_PNO_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_PNO_CCSD_PNO_H_

#include <iostream>
#include <sstream>
#include <tiledarray.h>

#include "mpqc/chemistry/qc/scf/mo_build.h"
#include "mpqc/chemistry/qc/mbpt/denom.h"
#include "mpqc/chemistry/qc/wfn/lcao_wfn.h"
#include "mpqc/math/tensor/clr/decomposed_tensor_algebra.h"

namespace mpqc {
namespace lcao {

  /**
   *  \brief CCSD_PNO class
   *
   *  keyword to call this class CCSD-PNO
   */

  template <typename Tile, typename Policy>
  class CCSD_PNO : public LCAOWavefunction<Tile, Policy> {

  private:
   //const KeyVal kv_;
   bool df_;
   double tcut_;
   std::shared_ptr<Wavefunction> ref_wfn_;
   TA::DistArray<Tile, Policy>  t2_mp2_;

   void init();

   // compute MP2 T2 amplitudes
   void compute_mp2_t2();
   // compute the occ and vir matrices for reblocking T2
   void compute_M_reblock(TA::DistArray<Tile, Policy> &occ_convert,
                          TA::DistArray<Tile, Policy> &vir_convert);
   // svd decompostion of T2
   void pno_decom();


  public:
   CCSD_PNO() = default;
   CCSD_PNO(const KeyVal &kv);

   ~CCSD_PNO() = default;
   void compute(PropertyBase* pb) override;
   // compute function
   double value() override;

  };

//#if TA_DEFAULT_POLICY == 0
//extern template class CCSD_PNO<TA::TensorD, TA::DensePolicy>;
//#elif TA_DEFAULT_POLICY == 1
//extern template class CCSD_PNO<TA::TensorD, TA::SparsePolicy>;
//#endif

}  // namespace lcao
}  // namespace mpqc

#include "ccsd_pno_impl.h"

#endif // MPQC4_SRC_MPQC_CHEMISTRY_QC_PNO_CCSD_PNO_H_
