//
// Created by Chong Peng on 3/2/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_FACTORY_UTILITY_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_FACTORY_UTILITY_H_

#include <cwchar>
#include <iostream>
#include <string>
#include <vector>

#include <libint2/engine.h>

#include "mpqc/chemistry/molecule/molecule.h"
#include "mpqc/chemistry/qc/lcao/basis/basis_registry.h"
#include "mpqc/chemistry/qc/lcao/expression/formula.h"
#include "mpqc/chemistry/qc/lcao/expression/orbital_registry.h"
#include "mpqc/chemistry/qc/lcao/integrals/integrals.h"
#include "mpqc/chemistry/qc/lcao/integrals/make_engine.h"
#include "mpqc/chemistry/qc/lcao/integrals/task_integrals.h"
#include "mpqc/util/meta/make_array.h"
#include "mpqc/util/misc/pool.h"

namespace mpqc {
namespace lcao {
namespace gaussian {

namespace detail {

libint2::Operator to_libint2_operator(Operator::Type mpqc_oper);

libint2::any to_libint2_operator_params(
    Operator::Type mpqc_oper, const Molecule &molecule,
    const std::vector<std::pair<double, double>> &gtg_params =
        std::vector<std::pair<double, double>>());

/// make screener
std::shared_ptr<Screener> make_screener(madness::World &world,
                                        const ShrPool<libint2::Engine> &engine,
                                        const BasisVector &bases,
                                        const std::string &screen,
                                        double screen_threshold);

/// given OrbitalIndex, find the corresponding basis
std::shared_ptr<Basis> index_to_basis(
    const OrbitalBasisRegistry &basis_registry, const OrbitalIndex &index);

///
/// Functions to generate formula
///

// clang-format off
/**
  * Given Formula with rank = 4, return DensityFitting formula
  *
  * This function is also used in LCAOFactory density fitting formula
  *parsing
  *
  * @param formula that has string format (in1 in2 | oper | in3 in4 ) or <in1 in2 | oper | in3 in4 >
  * @return array of wstring with length 3
  *         - string0 ( dfbs | oper | in1 in2 )[inv_sqr]  or ( dfbs |oper | in1 in3 ) [inv_sqr]
  *         - string1 ( dfbs | oper | in3 in4 )[inv_sqr]  or ( dfbs | oper | in2 in4 )[inv_sqr]
  *
  *         where ( dfbs | oper | in1 in2 )[inv_sqr]  = ( dfbs | oper | dfbs )[inv_sqr] * ( dfbs | oper | in1 in2 )
  *
  **/
// clang-format on
std::array<std::wstring, 2> get_df_formula(const Formula &formula);

// clang-format off
/**
   *  Given formula with rank = 2 and J or K operation, return the G integral
   *
   *  @param Formula that has string format ( in1 | oper | in2 ) or < in1 | oper| in2 > where oper is J or K operation
   *  @return Formula
   *          - Formula that has string format ( in1 in2 | G | obs obs ) or <in1 obs| G | in2 obs > for J
   *          - Formula that has string format ( in1 obs | G | in2 obs ) or <in1 in2| G | obs obs > for K, KAlpha, KBeta
   */
// clang-format on
Formula get_jk_formula(const Formula &formula, const std::wstring &obs);

// clang-format off
/**
  * Given formula with rank = 2 and J or K operation, return the G integral
  *with DensityFitting
  *
  * @param Formula that has string format ( in1 | oper | in2 ) or < in1 | oper | in2 >, where oper is J or K operation
  * @return result array of Formula with size 3
  *         - 3 Formula that has string format ( dfbs | G | in1 in2 ) ( dfbs | G | dfbs )[inv] ( dfbs | G | obs obs ) for J
  *         - 3 Formula that has string format ( dfbs | G | in1 obs ) ( dfbs | G | dfbs )[inv] ( dfbs | G | in2 obs ) for K, KAlpha, KBeta
  */
// clang-format on
std::array<Formula, 3> get_jk_df_formula(const Formula &formula,
                                         const std::wstring &obs);

// clang-format off
/**
  * Given formula with rank = 2 and Fock operation, return 3 formula to compute
  *it
  *
  * @param Formula that has string format ( in1 | oper | in2 ), where oper is
  *Fock operation
  * @return array of Formula with size 3
  *         - 3 Formula that has string format ( in1 | H | in2 ) ( in1 | J | in2 ) ( in1 |K | in2 ) for Fock
  *         - 3 Formula that has string format ( in1 | H | in2 ) ( in1 | J | in2 ) ( in1 | KAlpha | in2 ) for FockAlpha
  *         - 3 Formula that has string format ( in1 | H | in2 ) ( in1 | J | in2 ) ( in1 | KBeta | in2 ) for FockBeta
  */
// clang-format on
std::array<Formula, 3> get_fock_formula(const Formula &formula);

// clang-format off
/**
 * Given operation that is J or K operation, return the orbital index that
 *maps to density
 * @param operation J or K operation
 *
 * @return result OrbitalIndex
 *         - m for J or K
 *         - m_α for KAlpha
 *         - m_β for KBeta
 */
// clang-format on
OrbitalIndex get_jk_orbital_space(const Operator &operation);

/*!
 * \brief This takes real or imaginary part from a complex array
 * \tparam Policy can be TA::SparsePolicy or TA::DensePolicy
 * \tparam Is_real true if user requests real part, false if imaginary
 * \param a TensorZ array
 * \return a TensorD array
 */
template <typename Policy, bool is_real = true>
TA::DistArray<TA::TensorD, Policy> tensorZ_to_tensorD(
    const TA::DistArray<TA::TensorZ, Policy> &complex_array) {
  TA::DistArray<TA::TensorD, Policy> result_array;

  auto take_part_from_tile = [=](TA::TensorD &result_tile,
                                 const TA::TensorZ &arg_tile) {
    const auto &range = arg_tile.range();
    const auto volume = range.volume();
    result_tile = TA::TensorD(range);

    float norm = 0.0;
    for (auto ord = 0; ord < volume; ++ord) {
      const auto z = arg_tile[ord];
      auto result_el = (is_real) ? z.real() : z.imag();
      norm += result_el * result_el;
      result_tile[ord] = result_el;
    }

    return std::sqrt(norm);
  };

  result_array =
      TA::foreach<TA::TensorD, TA::TensorZ>(complex_array, take_part_from_tile);

  complex_array.world().gop.fence();

  return result_array;
}

/*!
 * \brief This computes non-negligible shell pair list; shells \c i and \c j
 * form a non-negligible pair if they share a center or the Frobenius norm of
 * their overlap is greater than threshold
 * \param basis1 a basis
 * \param basis2 a basis
 * \param threshold
 *
 * \return a list of pairs with
 * key: shell index
 * mapped value: a vector of shell indices
 */
std::vector<std::vector<size_t>> parallel_compute_shellpair_list(
    madness::World &world, const Basis &basis1, const Basis &basis2,
    double threshold = 1e-12, double engine_precision = 0.0);

}  // namespace detail

}  // namespace gaussian

namespace detail {
/// find the corresponding AO formula, if index is already AO, it will be
/// ignored
template <typename Array>
Formula lcao_to_ao(const Formula &formula,
                   const OrbitalSpaceRegistry<Array> &orbital_registry) {
  std::vector<OrbitalIndex> ao_left_index, ao_right_index;

  int increment = 0;
  auto left_index = formula.bra_indices();
  for (const auto &index : left_index) {
    // find the correspoding ao index
    if (index.is_lcao()) {
      auto ao_index = orbital_registry.retrieve(index).ao_index().name();
      ao_index = ao_index + std::to_wstring(increment);
      ao_left_index.push_back(ao_index);
      increment++;
    }
    // if already ao, do nothing
    else {
      ao_left_index.push_back(index);
    }
  }

  auto right_index = formula.ket_indices();
  for (const auto &index : right_index) {
    // find the correspoding ao index
    if (index.is_lcao()) {
      auto ao_index = orbital_registry.retrieve(index).ao_index().name();
      ao_index = ao_index + std::to_wstring(increment);
      ao_right_index.push_back(ao_index);
      increment++;
    }
    // if already ao, do nothing
    else {
      ao_right_index.push_back(index);
    }
  }

  // set formula with ao index
  auto ao_formula = formula;
  ao_formula.set_bra_indices(ao_left_index);
  ao_formula.set_ket_indices(ao_right_index);

  return ao_formula;
}

/// check if all index in formula are in LCAO
bool if_all_lcao(const Formula &formula);

/// check if all index in formula are in AO
bool if_all_ao(const Formula &formula);
}

}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_FACTORY_UTILITY_H_
