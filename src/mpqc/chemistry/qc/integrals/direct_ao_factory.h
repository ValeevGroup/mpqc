//
// Created by Chong Peng on 9/12/16.
//

#ifndef MPQC_DIRECT_AO_FACTORY_H
#define MPQC_DIRECT_AO_FACTORY_H

#include <type_traits>

#include "mpqc/util/misc/time.h"
#include <madness/world/worldmem.h>
#include "mpqc/chemistry/qc/integrals/ao_factory_base.h"

namespace mpqc {
namespace integrals {

template <typename Tile, typename Policy>
class DirectAOFactory;

namespace detail{
template <typename Tile, typename Policy>
std::shared_ptr<DirectAOFactory<Tile,Policy>> construct_direct_ao_factory(const KeyVal& kv){
  std::shared_ptr<DirectAOFactory<Tile,Policy>> direct_ao_factory;
  if(kv.exists_class("wfn_world:direct_ao_factory")){
    direct_ao_factory = kv.class_ptr<DirectAOFactory<Tile,Policy>>("wfn_world:direct_ao_factory");
  }
  else{
    direct_ao_factory = std::make_shared<DirectAOFactory<Tile,Policy>>(kv);
    std::shared_ptr<DescribedClass> direct_ao_factory_base = direct_ao_factory;
    KeyVal& kv_nonconst = const_cast<KeyVal&>(kv);
    kv_nonconst.keyval("wfn_world").assign("direct_ao_factory",direct_ao_factory_base);
  }
  return direct_ao_factory;
};
}


/**
 * \brief Direct AOFactory Class
 *
 *
 */


template <typename Tile, typename Policy>
class DirectAOFactory : public AOFactoryBase, public DescribedClass {
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
  DirectAOFactory() = default;
  DirectAOFactory(DirectAOFactory&&) = default;

  DirectAOFactory& operator=(DirectAOFactory&&) = default;

  /**
   * KeyVal constructor
   * It takes all the options from AOFactoryBase
   * @param accurate_time, bool, control if use fence in timing, default false
   *
   * @return
   */

  DirectAOFactory(const KeyVal& kv)
      : AOFactoryBase(kv), direct_ao_formula_registry_() {

    std::string prefix = "";
    if(kv.exists("wfn_wolrd") || kv.exists_class("wfn_world")){
      prefix = "wfn_world:";
    }

    accurate_time_ = kv.value(prefix + "accurate_time", false);

    /// For other Tile type, need to implement set_oper();
    set_oper(Tile());
  }

  virtual ~DirectAOFactory() noexcept = default;

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
typename DirectAOFactory<Tile, Policy>::DirectTArray
DirectAOFactory<Tile, Policy>::compute(const Formula& formula) {
  DirectTArray result;

  // find if in registry
  auto iter = direct_ao_formula_registry_.find(formula);

  if (iter != direct_ao_formula_registry_.end()) {
    result = *(iter->second);
    utility::print_par(world_, "Retrieved Direct AO Integral: ");
    utility::print_par(world_, utility::to_string(formula.string()));
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
typename DirectAOFactory<Tile, Policy>::DirectTArray
DirectAOFactory<Tile, Policy>::compute2(const Formula& formula) {
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
    utility::print_par(world_, utility::to_string(formula.string()));
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
    utility::print_par(world_, utility::to_string(formula.string()));
    double size = mpqc::detail::array_size(result.array());
    utility::print_par(world_, " Size: ", size, " GB");
    utility::print_par(world_, " Time: ", time, " s\n");
  } else {
    throw std::runtime_error(
        "Unsupported Operator in DirectAOFactory!!\n");
  }

  madness::print_meminfo(
      world_.rank(),
      utility::wconcat("DirectAOFactory:", formula.string()));
  return result;
}

template <typename Tile, typename Policy>
typename DirectAOFactory<Tile, Policy>::DirectTArray
DirectAOFactory<Tile, Policy>::compute3(const Formula& formula) {
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
  utility::print_par(world_, utility::to_string(formula.string()));
  double size = mpqc::detail::array_size(result.array());
  utility::print_par(world_, " Size: ", size, " GB");
  utility::print_par(world_, " Time: ", time, " s\n");
  madness::print_meminfo(
      world_.rank(),
      utility::wconcat("DirectAOFactory:", formula.string()));

  return result;
}

template <typename Tile, typename Policy>
typename DirectAOFactory<Tile, Policy>::DirectTArray
DirectAOFactory<Tile, Policy>::compute4(const Formula& formula) {

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
  utility::print_par(world_, utility::to_string(formula.string()));
  double size = mpqc::detail::array_size(result.array());
  utility::print_par(world_, " Size: ", size, " GB");
  utility::print_par(world_, " Time: ", time, " s\n");
  madness::print_meminfo(world_.rank(),
                         utility::wconcat("AOFactory:", formula.string()));
  return result;
}

extern template
class DirectAOFactory<TA::TensorD, TA::SparsePolicy>;
//extern template
//class DirectAOFactory<TA::TensorD, TA::DensePolicy>;

}  // end of namespace integrals
}  // end of namespace mpqc

#endif  // MPQC_DIRECT_AO_FACTORY_H
