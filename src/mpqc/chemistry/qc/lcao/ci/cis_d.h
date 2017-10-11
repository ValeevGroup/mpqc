//
// Created by Chong Peng on 10/2/17.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_LCAO_CI_CIS_D_H_
#define SRC_MPQC_CHEMISTRY_QC_LCAO_CI_CIS_D_H_

#include "mpqc/chemistry/qc/lcao/ci/cis.h"
#include "mpqc/chemistry/qc/lcao/mbpt/denom.h"

namespace mpqc {
namespace lcao {

/**
 *   CIS(D) Class, to be implemented
 *
 *   Helmich, B.; Hättig, C. J. Chem. Phys. 2011, 135 (21), 214106.
 */

template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cis_d_double_amplitude(
    LCAOFactoryBase<Tile, Policy>& lcao_factory,
    const TA::DistArray<Tile, Policy>& cis_ampl,
    typename Tile::numeric_type cis_energy, bool df) {
  TA::DistArray<Tile, Policy> result;
  TA::DistArray<Tile, Policy> tmp;
  std::wstring postfix = df ? L"[df]" : L"";

  auto g_ijka = lcao_factory.compute(L"<i j|G|k a>" + postfix);

  tmp("a,b,i,j") = g_ijka("i,j,k,b") * cis_ampl("k,a");
  result("a,b,i,j") = -(tmp("a,b,i,j") + tmp("b,a,j,i"));

  if (df) {
    auto Xia = lcao_factory.compute(L"( Λ |G|i a)[inv_sqr]");
    auto Xab = lcao_factory.compute(L"( Λ |G|a b)[inv_sqr]");

    tmp("a,b,i,j") = Xab("K,c,b") * cis_ampl("j,c") * Xia("K,i,a");
    result("a,b,i,j") += tmp("a,b,i,j") + tmp("b,a,j,i");

  } else {
    auto g_iabc = lcao_factory.compute(L"<i a|G|b c>");

    TA::DistArray<Tile, Policy> tmp;
    tmp("a,b,i,j") = g_iabc("i,c,a,b") * cis_ampl("j,c");
    result("a,b,i,j") += tmp("a,b,i,j") + tmp("b,a,j,i");
  }

  auto f_ab = lcao_factory.compute(L"<a|F|b>" + postfix);
  auto f_ij = lcao_factory.compute(L"<i|F|j>" + postfix);

  EigenVector<typename Tile::numeric_type> eps_o =
      array_ops::array_to_eigen(f_ij).diagonal();
  EigenVector<typename Tile::numeric_type> eps_v =
      array_ops::array_to_eigen(f_ab).diagonal();

  EigenVector<typename Tile::numeric_type> ens(eps_o.rows() + eps_v.rows());
  ens << eps_o, eps_v;

  detail::d_abij_inplace(result, ens, eps_o.rows(),0,cis_energy);

  return result;
};

/**
 *  CIS(D) Class
 */

// template <typename Tile, typename Policy>
// class CIS_D : public CIS<Tile,Policy> {
//
//};

}  // namespace lcao
}  // namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_LCAO_CI_CIS_D_H_
