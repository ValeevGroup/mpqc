//
// Created by Chong Peng on 9/12/16.
//

#ifndef MPQC_DIRECT_ATOMIC_INTEGRAL_H
#define MPQC_DIRECT_ATOMIC_INTEGRAL_H

#include <type_traits>

#include "../../../../../utility/time.h"
#include <madness/world/worldmem.h>
#include <mpqc/chemistry/qc/integrals/atomic_integral_base.h>
#include <rapidjson/document.h>

namespace mpqc {
namespace integrals {

template <typename Tile, typename Policy>
class DirectAtomicIntegral;

namespace detail{
template <typename Tile, typename Policy>
std::shared_ptr<DirectAtomicIntegral<Tile,Policy>> construct_direct_atomic_integral(const KeyVal& kv){
  std::shared_ptr<DirectAtomicIntegral<Tile,Policy>> direct_ao_int;
  if(kv.exists_class("wfn_world:direct_atomic_integral")){
    direct_ao_int = kv.class_ptr<DirectAtomicIntegral<Tile,Policy>>("wfn_world:direct_atomic_integral");
  }
  else{
    direct_ao_int = std::make_shared<DirectAtomicIntegral<Tile,Policy>>(kv);
    std::shared_ptr<DescribedClass> direct_ao_int_base = direct_ao_int;
    KeyVal& kv_nonconst = const_cast<KeyVal&>(kv);
    kv_nonconst.keyval("wfn_world").assign("direct_atomic_integral",direct_ao_int_base);
  }
  return direct_ao_int;
};
}


/**
 * \brief Direct Atomic Integral Class
 *
 *
 */


template <typename Tile, typename Policy>
class DirectAtomicIntegral : public AtomicIntegralBase, public DescribedClass {
 public:
  using DirectTArray = integrals::DirectArray<Tile, Policy>;
  using TArray = TA::DistArray<integrals::DirectTile<Tile>, Policy>;
  using Op = std::function<Tile(TA::TensorD&&)>;

  /// private class member
 private:
  FormulaRegistry<DirectTArray> direct_ao_formula_registry_;
  Op op_;
  bool accurate_time_;

 public:
  /**
   * Default Constructor
   */
  DirectAtomicIntegral() = default;
  DirectAtomicIntegral(DirectAtomicIntegral&&) = default;

  DirectAtomicIntegral& operator=(DirectAtomicIntegral&&) = default;

  /**
   * Constructor
   *  @param world madness::World object
   *  @param op op is a function that will take TA::TensorD as argument and
   * return Tile
   *  @param mol shared pointer to Molecule
   *  @param obs shared pointer to OrbitalBasisRegistry
   *  @param gtg_params  parameters used in computing f12 integrals
   *  @param in rapidjson Document object
   *
   */
  DirectAtomicIntegral(
      madness::World& world, Op op,
      const std::shared_ptr<Molecule>& mol,
      const std::shared_ptr<basis::OrbitalBasisRegistry>& obs,
      const std::vector<std::pair<double, double>>& gtg_params =
          std::vector<std::pair<double, double>>(),
      const rapidjson::Document& in = rapidjson::Document())
      : AtomicIntegralBase(world, mol, obs, gtg_params, in),
        direct_ao_formula_registry_(),
        op_(op) {
    if (in.IsObject()) {
      accurate_time_ =
          in.HasMember("AccurateTime") ? in["AccurateTime"].GetBool() : false;
    } else {
      accurate_time_ = false;
    }
    utility::print_par(world, "AccurateTime: ", accurate_time_, "\n");
  }

  /**
   * KeyVal constructor
   * It takes all the options from AtomicIntegralBase
   * @param accurate_time, bool, control if use fence in timing, default false
   *
   * @return
   */

  DirectAtomicIntegral(const KeyVal& kv)
      : AtomicIntegralBase(kv), direct_ao_formula_registry_() {
    accurate_time_ = kv.value("accurate_time", false);

    /// For other Tile type, need to implement set_oper();
    set_oper(Tile());
  }

  virtual ~DirectAtomicIntegral() noexcept = default;

  /// set oper based on Tile type
  template<typename T = Tile>
  void set_oper(typename std::enable_if<std::is_same<T,TA::TensorD>::value, T>::type && t){
    op_ = TA::Noop<TA::TensorD,true>();
  }
  /// wrapper to compute function
  DirectTArray compute(const std::wstring& str) {
    auto formula = Formula(str);
    return compute(formula);
  }

  /**
   * compute direct integral by Formula
   */

  DirectTArray compute(const Formula&);

  /**
   * compute direct integral and return TA::expression
   */
  TA::expressions::TsrExpr<TArray, true> operator()(const std::wstring& str) {
    auto formula = Formula(str);
    auto array = compute(formula);
    auto& result = direct_ao_formula_registry_.retrieve(formula);
    return result.array()(formula.to_ta_expression());
  }

  /**
   * return direct ao formula registry
   */
  const FormulaRegistry<DirectTArray>& registry() const {
    return direct_ao_formula_registry_;
  }

  /**
   * return direct ao formula registry
   */
  FormulaRegistry<DirectTArray>& registry() {
    return direct_ao_formula_registry_;
  }

 protected:
  /// compute integrals that has two dimension
  DirectTArray compute2(const Formula& formula_string);

  /// compute integrals that has three dimension
  DirectTArray compute3(const Formula& formula_string);

  /// compute integrals that has two dimension
  DirectTArray compute4(const Formula& formula_string);

 private:
  /// compute sparse array
  template <typename U = Policy>
  DirectArray<Tile,
              typename std::enable_if<std::is_same<U, TA::SparsePolicy>::value,
                                      TA::SparsePolicy>::type>
  compute_integrals(
      madness::World& world, ShrPool<libint2::Engine>& engine,
      Bvector const& bases,
      std::shared_ptr<Screener> p_screen =
          std::make_shared<integrals::Screener>(integrals::Screener{})) {
    auto result = mpqc::integrals::direct_sparse_integrals(world, engine, bases,
                                                           p_screen, op_);
    return result;
  }

  /// compute dense array
  template <typename U = Policy>
  DirectArray<Tile,
              typename std::enable_if<std::is_same<U, TA::DensePolicy>::value,
                                      TA::DensePolicy>::type>
  compute_integrals(
      madness::World& world, ShrPool<libint2::Engine>& engine,
      Bvector const& bases,
      std::shared_ptr<Screener> p_screen =
          std::make_shared<integrals::Screener>(integrals::Screener{})) {
    auto result = mpqc::integrals::direct_dense_integrals(world, engine, bases,
                                                          p_screen, op_);
    return result;
  }
};

template <typename Tile, typename Policy>
typename DirectAtomicIntegral<Tile, Policy>::DirectTArray
DirectAtomicIntegral<Tile, Policy>::compute(const Formula& formula) {
  DirectTArray result;

  // find if in registry
  auto iter = direct_ao_formula_registry_.find(formula);

  if (iter != direct_ao_formula_registry_.end()) {
    result = *(iter->second);
    utility::print_par(world_, "Retrieved Direct AO Integral: ");
    utility::wprint_par(world_, formula.string());
    double size = mpqc::detail::array_size(result.array());
    utility::print_par(world_, " Size: ", size, " GB\n");
    return result;

  }else{
    // compute formula
    if (formula.rank() == 2) {
      result = compute2(formula);
      direct_ao_formula_registry_.insert(formula, result);
    } else if (formula.rank() == 3) {
      result = compute3(formula);
      direct_ao_formula_registry_.insert(formula, result);
    } else if (formula.rank() == 4) {
      result = compute4(formula);
      direct_ao_formula_registry_.insert(formula, result);
    }
    return result;
  }
}

template <typename Tile, typename Policy>
typename DirectAtomicIntegral<Tile, Policy>::DirectTArray
DirectAtomicIntegral<Tile, Policy>::compute2(const Formula& formula) {
  Bvector bs_array;
  double time = 0.0;
  std::shared_ptr<EnginePool<libint2::Engine>> engine_pool;
  mpqc::time_point time0;
  mpqc::time_point time1;

  DirectTArray result;

  if (formula.oper().is_onebody()) {
    time0 = mpqc::now(world_, accurate_time_);

    parse_one_body(formula, engine_pool, bs_array);
    result = compute_integrals(world_, engine_pool, bs_array);

    time1 = mpqc::now(world_, accurate_time_);
    time += mpqc::duration_in_s(time0, time1);

    utility::print_par(world_, "Computed Direct One Body Integral: ");
    utility::wprint_par(world_, formula.string());
    double size = mpqc::detail::array_size(result.array());
    utility::print_par(world_, " Size: ", size, " GB");
    utility::print_par(world_, " Time: ", time, " s\n");

  } else if (formula.oper().is_twobody()) {
    time0 = mpqc::now(world_, accurate_time_);

    parse_two_body_two_center(formula, engine_pool, bs_array);
    result = compute_integrals(world_, engine_pool, bs_array);

    time1 = mpqc::now(world_, accurate_time_);
    time += mpqc::duration_in_s(time0, time1);

    utility::print_par(world_, "Computed Direct Twobody Two Center Integral: ");
    utility::wprint_par(world_, formula.string());
    double size = mpqc::detail::array_size(result.array());
    utility::print_par(world_, " Size: ", size, " GB");
    utility::print_par(world_, " Time: ", time, " s\n");
  } else {
    throw std::runtime_error(
        "Unsupported Operator in DirectAtomicIntegral!!\n");
  }

  madness::print_meminfo(
      world_.rank(),
      utility::wconcat("DirectAtomicIntegral:", formula.string()));
  return result;
}

template <typename Tile, typename Policy>
typename DirectAtomicIntegral<Tile, Policy>::DirectTArray
DirectAtomicIntegral<Tile, Policy>::compute3(const Formula& formula) {
  double time = 0.0;
  mpqc::time_point time0;
  mpqc::time_point time1;
  time0 = mpqc::now(world_, accurate_time_);
  DirectTArray result;

  Bvector bs_array;
  std::shared_ptr<EnginePool<libint2::Engine>> engine_pool;
  std::shared_ptr<Screener> p_screener =
      std::make_shared<integrals::Screener>(integrals::Screener{});

  parse_two_body_three_center(formula, engine_pool, bs_array, p_screener);
  result = compute_integrals(this->world_, engine_pool, bs_array, p_screener);

  time1 = mpqc::now(world_, accurate_time_);
  time += mpqc::duration_in_s(time0, time1);

  utility::print_par(world_, "Computed Direct Twobody Three Center Integral: ");
  utility::wprint_par(world_, formula.string());
  double size = mpqc::detail::array_size(result.array());
  utility::print_par(world_, " Size: ", size, " GB");
  utility::print_par(world_, " Time: ", time, " s\n");
  madness::print_meminfo(
      world_.rank(),
      utility::wconcat("DirectAtomicIntegral:", formula.string()));

  return result;
}

template <typename Tile, typename Policy>
typename DirectAtomicIntegral<Tile, Policy>::DirectTArray
DirectAtomicIntegral<Tile, Policy>::compute4(const Formula& formula) {

  if(formula.notation() != Formula::Notation::Chemical){
    throw std::runtime_error("Direct AO Integral Only Support Chemical Notation! \n");
  }

  double time = 0.0;
  mpqc::time_point time0;
  mpqc::time_point time1;
  DirectTArray result;

  time0 = mpqc::now(world_, accurate_time_);

  Bvector bs_array;
  std::shared_ptr<Screener> p_screener =
      std::make_shared<integrals::Screener>(integrals::Screener{});
  std::shared_ptr<EnginePool<libint2::Engine>> engine_pool;

  parse_two_body_four_center(formula, engine_pool, bs_array, p_screener);

  result = compute_integrals(this->world_, engine_pool, bs_array, p_screener);

  time1 = mpqc::now(world_, accurate_time_);
  time += mpqc::duration_in_s(time0, time1);

  utility::print_par(world_, "Computed Direct Twobody Four Center Integral: ");
  utility::wprint_par(world_, formula.string());
  double size = mpqc::detail::array_size(result.array());
  utility::print_par(world_, " Size: ", size, " GB");
  utility::print_par(world_, " Time: ", time, " s\n");
  madness::print_meminfo(world_.rank(),
                         utility::wconcat("AtomicIntegral:", formula.string()));
  return result;
}

}  // end of namespace integrals
}  // end of namespace mpqc

#endif  // MPQC_DIRECT_ATOMIC_INTEGRAL_H
