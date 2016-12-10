//
// Created by Chong Peng on 3/2/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_AO_FACTORY_BASE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_AO_FACTORY_BASE_H_

#include <cwchar>
#include <iostream>
#include <string>
#include <vector>

#include "integrals.h"
#include "mpqc/chemistry/molecule/molecule.h"
#include "mpqc/chemistry/qc/basis/basis_registry.h"
#include "mpqc/chemistry/qc/expression/formula.h"
#include "mpqc/chemistry/qc/expression/orbital_registry.h"
#include "mpqc/chemistry/qc/integrals/make_engine.h"
#include "mpqc/chemistry/qc/integrals/task_integrals.h"
#include "mpqc/util/meta/make_array.h"

#include <libint2/engine.h>
#include "mpqc/util/misc/pool.h"

namespace mpqc {
namespace lcao {
namespace gaussian {

/**
 *
 * \brief base class for AOFactory
 *
 *  Options in Input
 *  @param Screen, string, name of screen method to use, default none
 *  @param Threshold, double, screen threshold, qqr or schwarz, default 1.0e-10
 *  @param Precision, double, precision in computing integral, default
 * std::numeric_limits<double>::epsilon()
   *
 */

class AOFactoryBase {
 public:
  using gtg_params_t = std::vector<std::pair<double, double>>;

  AOFactoryBase() = default;

  // clang-format off
  /**
   * \brief  KeyVal constructor
   *
   * It takes all the keys to construct OrbitalBasisRegistry and also the following
   *
   *  | KeyWord | Type | Default| Description |
   *  |---------|------|--------|-------------|
   *  |molecule|Molecule|none|keyval to construct molecule|
   *  |screen|string|none|method of screening, qqr or schwarz |
   *  |threshold|double|1e-10| screening threshold |
   *  |precision|double|std::numeric_limits<double>::epsilon() | integral precision |
   *  |corr_functions|int|6|f12 n of corr function,valid if aux_basis exsist in OrbitalBasisRegistry|
   *  |corr_param|int|0|f12 corr param, ,valid if aux_basis exsist in OrbitalBasisRegistry|
   */
  // clang-format on
  AOFactoryBase(const KeyVal &kv);

  virtual ~AOFactoryBase() = default;

  /// @return MADNESS world
  madness::World &world() { return world_; }

  /// @brief Molecule accessor
  /// @return molecule object
  const Molecule &molecule() const { return *mol_; }

  /// @brief (contracted) Gaussian-types geminal parameters accessor
  /// @return Gaussian-type geminal parameters
  const gtg_params_t &gtg_params() const {
    TA_USER_ASSERT(not gtg_params_.empty(),
                   "Gaussian-type geminal not initialized");
    return gtg_params_;
  }

  /// set OrbitalBasisRegistry
  void set_orbital_basis_registry(
      const std::shared_ptr<gaussian::OrbitalBasisRegistry> &obs) {
    orbital_basis_registry_ = obs;
  }

  /// @return the OrbitalBasisRegistry object
  const OrbitalBasisRegistry &orbital_basis_registry() const {
    return *orbital_basis_registry_;
  }

  OrbitalBasisRegistry &orbital_basis_registry() {
    return *orbital_basis_registry_;
  }

  const std::shared_ptr<OrbitalBasisRegistry>
  orbital_basis_registry_ptr() const {
    return orbital_basis_registry_;
  }
  /**
    * Given Formula with rank = 4, return DensityFitting formula
    *
    * This function is also used in LCAOFactory density fitting formula
    *parsing
    *
    * @param formula that has string format (in1 in2 | oper | in3 in4 ) or <in1
    *in2 | oper | in3 in4 >
    * @return array of wstring with length 3
    *         - string0 ( dfbs | oper | in1 in2 )  or ( dfbs |oper | in1 in3 )
    *         - string1 ( dfbs | oper | dfbs )[inv]
    *         - string2 ( dfbs | oper | in3 in4 )  or ( dfbs | oper | in2 in4 )
    */
  std::array<std::wstring, 3> get_df_formula(const Formula &formula);

  /// given OrbitalIndex, find the correspoding basis
  std::shared_ptr<Basis> index_to_basis(const OrbitalIndex &index) {
    auto iter = orbital_basis_registry_->find(index);
    if (iter == orbital_basis_registry_->end()) {
      throw std::runtime_error("Basis Set Not Found!!");
    } else {
      return iter->second;
    }
  }

 protected:
  /// parse operation and return one body engine
  libint2::Engine make_engine(const Operator &oper, int64_t max_nprim,
                              int64_t max_am);

  /// parse one body formula and set engine_pool and basis array
  void parse_one_body(const Formula &formula,
                      std::shared_ptr<utility::TSPool<libint2::Engine>> &engine_pool,
                      BasisVector &bases);

  /// parse two body two center formula and set two body kernel, basis array,
  /// max_nprim and max_am
  void parse_two_body_two_center(
      const Formula &formula,
      std::shared_ptr<utility::TSPool<libint2::Engine>> &engine_pool,
      BasisVector &bases);

  /// parse two body three center formula and set two body kernel, basis array,
  /// max_nprim and max_am
  void parse_two_body_three_center(
      const Formula &formula,
      std::shared_ptr<utility::TSPool<libint2::Engine>> &engine_pool, BasisVector &bases,
      std::shared_ptr<Screener> &p_screener);

  /// parse two body four center formula and set two body kernel, basis array,
  /// max_nprim and max_am
  void parse_two_body_four_center(
      const Formula &formula,
      std::shared_ptr<utility::TSPool<libint2::Engine>> &engine_pool, BasisVector &bases,
      std::shared_ptr<Screener> &p_screener);

  std::shared_ptr<Screener> make_screener_four_center(
      ShrPool<libint2::Engine> &engine, Basis &basis);

  std::shared_ptr<Screener> make_screener_three_center(
      ShrPool<libint2::Engine> &engine, Basis &dfbasis,
      Basis &basis);

  /**
   *  Given formula with rank = 2 and J or K operation, return the G integral
   *
   *  @param Formula that has string format ( in1 | oper | in2 ) or < in1 | oper
   *| in2 > where oper is J or K operation
   *  @return Formula
   *          - Formula that has string format ( in1 in2 | G | obs obs ) or <in1
   *obs| G | in2 obs > for J
   *          - Formula that has string format ( in1 obs | G | in2 obs ) or <in1
   *in2| G | obs obs > for K, KAlpha, KBeta
   */
  Formula get_jk_formula(const Formula &formula, const std::wstring &obs);

  /**
   * Given formula with rank = 2 and J or K operation, return the G integral
   *with DensityFitting
   *
   * @param Formula that has string format ( in1 | oper | in2 ) or < in1 | oper
   *| in2 >, where oper is J or K operation
   * @return result array of Formula with size 3
   *         - 3 Formula that has string format ( dfbs | G | in1 in2 ) ( dfbs |
   *G | dfbs )[inv] ( dfbs | G | obs obs ) for J
   *         - 3 Formula that has string format ( dfbs | G | in1 obs ) ( dfbs |
   *G | dfbs )[inv] ( dfbs | G | in2 obs ) for K, KAlpha, KBeta
   */
  std::array<Formula, 3> get_jk_df_formula(const Formula &formula,
                                           const std::wstring &obs);

  /**
   * Given formula with rank = 2 and Fock operation, return 3 formula to compute
   *it
   *
   * @param Formula that has string format ( in1 | oper | in2 ), where oper is
   *Fock operation
   * @return array of Formula with size 3
   *         - 3 Formula that has string format ( in1 | H | in2 ) ( in1 | J |
   *in2 ) ( in1 |K | in2 ) for Fock
   *         - 3 Formula that has string format ( in1 | H | in2 ) ( in1 | J |
   *in2 ) ( in1 | KAlpha | in2 ) for FockAlpha
   *         - 3 Formula that has string format ( in1 | H | in2 ) ( in1 | J |
   *in2 ) ( in1 | KBeta | in2 ) for FockBeta
   */
  std::array<Formula, 3> get_fock_formula(const Formula &formula);

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
  OrbitalIndex get_jk_orbital_space(const Operator &operation);

 protected:
  madness::World &world_;
  std::shared_ptr<OrbitalBasisRegistry> orbital_basis_registry_;

  // TODO these specify operator params, need to abstract out better
  std::shared_ptr<Molecule> mol_;
  gtg_params_t gtg_params_;
  std::string screen_;
  double screen_threshold_;
  double precision_;
};

namespace detail {

libint2::Operator to_libint2_operator(Operator::Type mpqc_oper);

libint2::any to_libint2_operator_params(Operator::Type mpqc_oper,
                                        const AOFactoryBase &base);
}

}  // namespace gaussian
}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_AO_FACTORY_BASE_H_
