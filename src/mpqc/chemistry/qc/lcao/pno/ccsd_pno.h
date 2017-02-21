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

   // compute MP2 T2 amplitudes
   TA::DistArray<Tile, Policy> compute_mp2_t2();
   // compute occ and vir matrices for reblocking MP2 T2
   void compute_M_reblock(TA::DistArray<Tile, Policy> &occ_convert,
                          TA::DistArray<Tile, Policy> &vir_convert);

   // compute CCSD with explicit initial values for T2
   // based on Chong's compute_ccsd_df(TArray &t1, TArray &t2) function
   double compute_ccsdpno_df(TA::DistArray<Tile, Policy> &t1,
                             TA::DistArray<Tile, Policy> &t2);

   // function for testing PNO decomposition ideas
   double pno_simul();
   // decompose T2 using PNO (eigen) decomposition or SVD
   void decom_t2(TA::DistArray<Tile, Policy> &t2_mp2, const DecomType = eigen_decom);
   // decompose T2 with MP2 and CCSD components
   void decom_t2(TA::DistArray<Tile, Policy> &t2_mp2,
                 const TA::DistArray<Tile, Policy> &t2_ccsd);

   // compute PNO coefficients: d^ij_a aij
   // where each tile is the size of a * aij (aij is PNO)
   TA::DistArray<Tile, Policy>
   compute_pno_coef(const TA::DistArray<Tile, Policy> &t2_mp2);
   // obtain the PNO coefficients in a vector
   // each vector entry contain ab matrix
   void comput_vec_pnocoef(const TA::DistArray<Tile, Policy>& t2_mp2,
                           const TA::TiledRange1& trange1_a);
   // construct vector<TArray> for ab^ij tensor (T2^ab_ij, g^ab_ij, and d^ab_ij)
   // where each entry contain reblocked ab matrix
   void comput_vecij_ab(const TA::DistArray<Tile, Policy>& ab_ij,
                         const TA::TiledRange1& trange1_a,
                         std::vector<TA::DistArray<Tile, Policy>>& vecij_ab,
                         const bool compute_pno_coef = true);
   // compute PNO transformed F^a_b
   void get_fab_pno(std::vector<TA::DistArray<Tile, Policy>>& vecij_fab);
   // compute PNO transformed G^ab_cd integrals
   void get_abcd_pno(std::vector<TA::DistArray<Tile, Policy>>& vecij_gabcd);
   // compute PNO transformed T2^ab_ij amplitudes or G^ab_ij integrals
   void get_abij_pno(const TA::DistArray<Tile, Policy>& ab_ij,
                     const TA::TiledRange1& trange1_a,
                     std::vector<TA::DistArray<Tile, Policy>>& vecij_pno);
   // compute CCD
   double compute_ccd_df();
   // compute PNO CCD
   double compute_ccd_pno();

  public:
   CCSD_PNO() = default;
   CCSD_PNO(const KeyVal &kv);

   ~CCSD_PNO() = default;

  protected:
   // compute function
   void evaluate(Energy* result) override;
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
