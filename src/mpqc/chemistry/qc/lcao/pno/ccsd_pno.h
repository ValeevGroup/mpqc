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

#include "mpqc/chemistry/qc/lcao/pno/util.h"
#include "mpqc/chemistry/qc/lcao/cc/ccsd.h"
#include "mpqc/chemistry/qc/lcao/mbpt/denom.h"
#include "mpqc/chemistry/qc/lcao/scf/mo_build.h"
#include "mpqc/chemistry/qc/lcao/wfn/lcao_wfn.h"
#include "mpqc/math/tensor/clr/decomposed_tensor_algebra.h"
#include "mpqc/chemistry/qc/lcao/expression/trange1_engine.h"

namespace mpqc {
namespace lcao {

  enum DecomType {eigen_decom, svd};

  /**
   *  \brief CCSD_PNO class
   *
   *  keyword to call this class CCSD-PNO
   */

  template <typename Tile, typename Policy>
  class CCSD_PNO : public CCSD<Tile, Policy> {

  private:
   double tcut_;
   std::vector<TA::DistArray<Tile, Policy>> vecij_dab_;

   /// compute MP1 T2 amplitudes
   TA::DistArray<Tile, Policy> compute_mp2_t2();
   /** compute occ and vir matrices for reblocking T2 amplitudes:
    * occ_convert has size of old_i * new_i
    * vir_convert is old_a * new_a
    */
   void compute_M_reblock(TA::DistArray<Tile, Policy> &occ_convert,
       TA::DistArray<Tile, Policy> &vir_convert);

   /// compute CCSD with explicit initial values for T2
   /// based on Chong's compute_ccsd_df(TArray &t1, TArray &t2) function
   double compute_ccsdpno_df(TA::DistArray<Tile, Policy> &t1,
       TA::DistArray<Tile, Policy> &t2);

   /** PNO simulation and test PNO decomposition ideas
    *  first, initial T2 amplitudes for CCSD iterations are decomposed
    *  second, those amplitudes are reconstructed with some truncations
    *  finally, the truncated amplitudes are used as initial guess for CCSD iterations.
    */
   double pno_simul();
   /// decompose T2 using PNO (eigen) decomposition or SVD
   void decom_t2(TA::DistArray<Tile, Policy> &t2_mp2, const DecomType = eigen_decom);

  public:
   CCSD_PNO() = default;

   /**
    * KeyVal constructor
    * @param kv
    *
    * keywords : all keywords for LCAOWavefunciton
    *
    * | KeyWord | Type | Default| Description |
    * |---------|------|--------|-------------|
    * | ref     | wfn  | none   | reference wavefunction, RHF for example |
    * | df      | bool | false  | choice of using density fitting
    */
   CCSD_PNO(const KeyVal &kv);

   ~CCSD_PNO() = default;

  protected:
   /// compute function
   void evaluate(Energy* result) override;
  };

#if TA_DEFAULT_POLICY == 0
extern template class CCSD_PNO<TA::TensorD, TA::DensePolicy>;
#elif TA_DEFAULT_POLICY == 1
extern template class CCSD_PNO<TA::TensorD, TA::SparsePolicy>;
#endif

}  // namespace lcao
}  // namespace mpqc

#include "ccsd_pno_impl.h"

#endif // MPQC4_SRC_MPQC_CHEMISTRY_QC_PNO_CCSD_PNO_H_
