//
// Created by Chong Peng on 1/7/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_LCAO_FACTORY_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_LCAO_FACTORY_H_

#include <string>
#include <vector>

#include <tiledarray.h>


#include "mpqc/math/linalg/diagonal_array.h"
#include "mpqc/chemistry/qc/lcao/expression/orbital_registry.h"
#include "mpqc/chemistry/qc/lcao/integrals/ao_factory.h"
#include "mpqc/chemistry/qc/lcao/wfn/wfn_world.h"


namespace mpqc {
namespace lcao {

template <typename Tile, typename Policy>
class LCAOFactory;


namespace detail{


template <typename Tile, typename Policy>
std::shared_ptr<LCAOFactory<Tile,Policy>> construct_lcao_factory(const KeyVal& kv){
  std::shared_ptr<LCAOFactory<Tile,Policy>> lcao_factory;
  if(kv.exists_class("wfn_world:lcao_factory")){
    lcao_factory = kv.class_ptr<LCAOFactory<Tile,Policy>>("wfn_world:lcao_factory");
  }
  else{
    lcao_factory = std::make_shared<LCAOFactory<Tile,Policy>>(kv);
    std::shared_ptr<DescribedClass> ao_factory_base = lcao_factory;
    KeyVal& kv_nonconst = const_cast<KeyVal&>(kv);
    kv_nonconst.keyval("wfn_world").assign("lcao_factory",ao_factory_base);
  }
  return lcao_factory;
};


}  // namespace  detail

// TODO MO transform that minimize operations by permutation
/**
 * \brief Molecule Integral computation class
 *  This class computes molecular integrals using a Formula object
 *
 *  compute(formula) return TArray object
 *  (formula) return TArray expression
 *
 */
template <typename Tile, typename Policy>
class LCAOFactory : virtual public DescribedClass{
 public:
  using TArray = TA::DistArray<Tile, Policy>;
  // for now hardwire to Gaussians
  // TODO generalize to non-gaussian AO operators
  using AOFactoryType = gaussian::AOFactory<Tile, Policy>;

  /**
   * Constructor
   * @param WfnWorld
   * @param KeyVal
   *
   * KeyVal options
   * @param accurate_time, if do fence when timing, default false
   * @param keep_partial_transform, if use strength reduction, default false
   *
   */
  LCAOFactory(const KeyVal& kv)
    : world_(*kv.value<madness::World *>("$:world")),
      ao_factory_(*gaussian::construct_ao_factory<Tile, Policy>(kv)),
      orbital_space_registry_(std::make_shared<OrbitalSpaceRegistry<TArray>>()),
      mo_formula_registry_()
  {
    std::string prefix = "";
    if(kv.exists("wfn_world") || kv.exists_class("wfn_world")){
      prefix = "wfn_world:";
    }
    accurate_time_ = kv.value<bool>(prefix + "accurate_time",false);
    keep_partial_transforms_ = kv.value<bool>(prefix + "keep_partial_transform",false);
    ao_factory_.set_orbital_space_registry(orbital_space_registry_);
  }

  void obsolete() {
    // obsolete self
    mo_formula_registry_.purge(world_);
    if(orbital_space_registry_!= nullptr){
      orbital_space_registry_->clear();
    }
    // obsolete AOFactory
    ao_factory_.obsolete();
  }

  /// return reference to madness::World
  madness::World& world() const { return world_; }

  /// return reference to AOFactory object
  AOFactoryType& ao_factory() const { return ao_factory_; }

  /// wrapper to operator() function in AOFactory
  TA::expressions::TsrExpr<TArray, true> ao_factory(
      const std::wstring& str) {
    return std::move(ao_factory_(str));
  };

  /// return OrbitalSpaceRegistry
  const OrbitalSpaceRegistry<TArray>& orbital_space() const {
    return *orbital_space_registry_;
  }

  OrbitalSpaceRegistry<TArray>& orbital_space() {
    return *orbital_space_registry_;
  }

  /// return reference to FormulaRegistry
  const FormulaRegistry<TArray>& registry() const {
    return mo_formula_registry_;
  }

  /// return reference to FormulaRegistry
  FormulaRegistry<TArray>& registry() { return mo_formula_registry_; }

  /// return accurate time
  bool accurate_time() const { return accurate_time_; }

  /// reports the partial tform flag; if true, partially-transformed integrals
  /// are stored
  /// @note at this time only supported for 3-index integrals
  bool keep_partial_transforms() const { return keep_partial_transforms_; }

  /// sets the partial tform flag; if true, partially-transformed integrals are
  /// stored
  /// @note at this time only supported for 3-index integrals
  void keep_partial_transforms(bool flag) { keep_partial_transforms_ = flag; }

  /// wrapper to compute function
  TArray compute(const std::wstring&);

  /**
   *  compute integral by Formula
   *  this function will look into registry first
   *  if Formula computed, it will return it from registry
   *  if not, it will compute it
   */
  TArray compute(const Formula&);

  /// compute with str and return expression
  TA::expressions::TsrExpr<TArray, true> operator()(const std::wstring& str) {
    auto formula = Formula(str);
    TArray array = compute(formula);
    auto& result = mo_formula_registry_.retrieve(formula);
    return result(formula.to_ta_expression());
  };

  /// purge formulae that contain Operator described by string \c str
  /// from mo_registry and ao_registry
  void purge_operator(madness::World& world, const std::wstring& str) {
    Operator oper(str);
    Operator::Type oper_type = oper.type();

    mo_formula_registry_.purge_operator(world, oper_type);
    ao_factory().registry().purge_operator(world, oper_type);
  }

  /// purge formulae that contain index described by string \c idx_str
  /// from mo_registry and ao_registry
  void purge_index(madness::World& world, const std::wstring& idx_str) {
    OrbitalIndex index(idx_str);
    mo_formula_registry_.purge_index(world, index);
    ao_factory().registry().purge_index(world, index);
  }

  /// purge formula described by string \c str
  /// from mo_registry
  void purge_formula(madness::World& world, const std::wstring& str) {
    mo_formula_registry_.purge_formula(world, str);
  }

 private:

  /// compute integrals that has two dimension
  TArray compute2(const Formula& formula_string);

  /// compute integrals that has three dimension
  TArray compute3(const Formula& formula_string);

  /// compute integrals that has four dimension
  TArray compute4(const Formula& formula_string);

 protected:
  /// find the corresponding AO formula, if index is already AO, it will be
  /// ignored
  Formula mo_to_ao(const Formula& formula);

  /// "reduce" by converting to AO at most 1 index in the formula,
  /// the index is chosen to minimize strength reduction (thus recursive
  /// reduction ensures maximum strength reduction)
  /// @return {"reduced" formula, {\c pos , \c rank }}, where \c pos and \c rank
  ///         are the location and rank of the reduced index
  std::tuple<Formula, std::pair<Formula::Position, size_t>> reduce_formula(
      const Formula& formula);

  /// assert all index in formula are in MO
  void assert_all_mo(const Formula& formula);

 protected:
  madness::World& world_;
  AOFactoryType& ao_factory_;
  std::shared_ptr<OrbitalSpaceRegistry<TArray>> orbital_space_registry_;
  FormulaRegistry<TArray> mo_formula_registry_;
  bool keep_partial_transforms_;  //!< if true, keep partially-transformed ints
                                  //!(false by default)
  bool accurate_time_;
};

template <typename Tile, typename Policy>
typename LCAOFactory<Tile, Policy>::TArray LCAOFactory<Tile, Policy>::compute2(
    const Formula& formula_string) {
  double time = 0.0;
  mpqc::time_point time0;
  mpqc::time_point time1;

  TArray result;
  // Identity matrix
  if (formula_string.oper().type() == Operator::Type::Identity) {
    time0 = mpqc::now(world_, accurate_time_);

    auto left_index1 = formula_string.bra_indices()[0];
    auto right_index1 = formula_string.ket_indices()[0];
    auto& left1 = orbital_space_registry_->retrieve(left_index1);
    auto& right1 = orbital_space_registry_->retrieve(right_index1);

    auto tr1 = left1.trange();
    auto tr2 = right1.trange();

    // create diagonal array
    result = array_ops::create_diagonal_array_from_eigen<Tile,Policy>(world_,tr1,tr2,1.0);
    result.truncate();

    time1 = mpqc::now(world_, accurate_time_);
    time += mpqc::duration_in_s(time0, time1);

    utility::print_par(world_, "Computed Identity: ",
                       utility::to_string(formula_string.string()));
    double size = mpqc::detail::array_size(result);
    utility::print_par(world_, " Size: ", size, " GB");
    utility::print_par(world_, " Time: ", time, " s\n");
    return result;
  }

  // get AO
  auto ao_formula = mo_to_ao(formula_string);
  auto ao_factory = ao_factory_.compute(ao_formula);

  time0 = mpqc::now(world_, accurate_time_);
  // convert to MO
  result = ao_factory;
  // get coefficient
  auto left_index1 = formula_string.bra_indices()[0];
  if (left_index1.is_mo()) {
    auto& left1 = orbital_space_registry_->retrieve(left_index1);
    result("i,r") = result("p,r") * left1("p,i");
  }
  auto right_index1 = formula_string.ket_indices()[0];
  if (right_index1.is_mo()) {
    auto& right1 = orbital_space_registry_->retrieve(right_index1);
    result("p,k") = result("p,r") * right1("r,k");
  }

  result.truncate();
  time1 = mpqc::now(world_, accurate_time_);
  time += mpqc::duration_in_s(time0, time1);
  utility::print_par(world_, "Transformed LCAO Integral: ",
                     utility::to_string(formula_string.string()));
  double size = mpqc::detail::array_size(result);
  utility::print_par(world_, " Size: ", size, " GB");
  utility::print_par(world_, " Time: ", time, " s\n");

  return result;
}

template <typename Tile, typename Policy>
typename LCAOFactory<Tile, Policy>::TArray LCAOFactory<Tile, Policy>::compute3(
    const Formula& formula_string) {
  double time = 0.0;
  mpqc::time_point time0;
  mpqc::time_point time1;

  TArray result;
  if (not keep_partial_transforms()) {  // compute from AO ints
    // get AO
    auto ao_formula = mo_to_ao(formula_string);
    auto ao_factory = ao_factory_.compute(ao_formula);

    time0 = mpqc::now(world_, accurate_time_);

    // transform to MO, only convert the right side
    // TODO optimize strength reduction,
    //      e.g. may need to transform last index first
    auto right_index1 = formula_string.ket_indices()[0];
    if (right_index1.is_mo()) {
      auto& right1 = orbital_space_registry_->retrieve(right_index1);
      result("K,i,q") = ao_factory("K,p,q") * right1("p,i");
    }
    auto right_index2 = formula_string.ket_indices()[1];
    if (right_index2.is_mo()) {
      auto& right2 = orbital_space_registry_->retrieve(right_index2);
      result("K,p,j") = result("K,p,q") * right2("q,j");
    }
  } else {  // tform to optimally reduce strength, store partial transform
            // results

    // compute reduced formula
    Formula reduced_formula;  // reduced formula
    std::pair<Formula::Position, size_t> reduced_index_coord;
    std::tie(reduced_formula, reduced_index_coord) =
        reduce_formula(formula_string);
    auto reduced_integral =
        reduced_formula.is_ao()
            ? ao_factory_.compute(reduced_formula)
            : (keep_partial_transforms() ? this->compute(reduced_formula)
                                         : this->compute3(reduced_formula));

    time0 = mpqc::now(world_, accurate_time_);

    const auto reduced_index_position = reduced_index_coord.first;
    const auto reduced_index_rank = reduced_index_coord.second;
    const auto reduced_index_absrank =
        reduced_index_rank + (reduced_index_position == Formula::Position::Bra
                                  ? 0
                                  : reduced_formula.bra_indices().size());

    // extract orbital coefficients for transformation
    auto reduced_index = (reduced_index_position == Formula::Position::Bra)
                             ? formula_string.bra_indices()[reduced_index_rank]
                             : formula_string.ket_indices()[reduced_index_rank];
    auto& reduced_index_space =
        orbital_space_registry_->retrieve(reduced_index);
    const auto& reduced_index_coeff = reduced_index_space.coefs();

    // transform
    auto result_key =
        formula_string.to_ta_expression(mpqc::detail::append_count(0));
    auto reduced_key =
        reduced_formula.to_ta_expression(mpqc::detail::append_count(0));
    auto coeff_key = reduced_index_space.ao_index().to_ta_expression() +
                     std::to_string(reduced_index_absrank) + ", " +
                     reduced_index.to_ta_expression() +
                     std::to_string(reduced_index_absrank);
    result(result_key) =
        reduced_integral(reduced_key) * reduced_index_coeff(coeff_key);
  }

  result.truncate();

  time1 = mpqc::now(world_, accurate_time_);
  time += mpqc::duration_in_s(time0, time1);

  utility::print_par(world_, "Transformed LCAO Integral: ",
                     utility::to_string(formula_string.string()));
  double size = mpqc::detail::array_size(result);
  utility::print_par(world_, " Size: ", size, " GB");
  utility::print_par(world_, " Time: ", time, " s\n");

  return result;
};

template <typename Tile, typename Policy>
typename LCAOFactory<Tile, Policy>::TArray LCAOFactory<Tile, Policy>::compute4(
    const Formula& formula_string) {
  double time = 0.0;
  mpqc::time_point time0;
  mpqc::time_point time1;
  TArray result;
  if (formula_string.oper().has_option(Operator::Option::DensityFitting)) {
    // get df formula
    auto df_formulas = ao_factory_.get_df_formula(formula_string);

    auto notation = formula_string.notation();
    // compute integral
    TArray left = compute(df_formulas[0]);

    TArray right = compute(df_formulas[2]);

    TArray center = ao_factory_.compute(df_formulas[1]);

    time0 = mpqc::now(world_, accurate_time_);

    if (notation == Formula::Notation::Chemical) {
      result("i,j,k,l") = left("q,i,j") * center("q,p") * right("p,k,l");
    } else {
      result("i,k,j,l") = left("q,i,j") * center("q,p") * right("p,k,l");
    }

    time1 = mpqc::now(world_, accurate_time_);
    time += mpqc::duration_in_s(time0, time1);

  } else {
    // get AO
    auto ao_formula = mo_to_ao(formula_string);
    auto ao_factory = ao_factory_.compute(ao_formula);

    // convert to MO
    time0 = mpqc::now(world_, accurate_time_);

    // get coefficient
    auto left_index1 = formula_string.bra_indices()[0];
    if (left_index1.is_mo()) {
      auto& left1 = orbital_space_registry_->retrieve(left_index1);
      result("i,q,r,s") = ao_factory("p,q,r,s") * left1("p,i");
    }


    auto left_index2 = formula_string.bra_indices()[1];
    if (left_index2.is_mo()) {
      auto& left2 = orbital_space_registry_->retrieve(left_index2);
      result("p,i,r,s") = result("p,q,r,s") * left2("q,i");
    }

    auto right_index1 = formula_string.ket_indices()[0];
    if (right_index1.is_mo()) {
      auto& right1 = orbital_space_registry_->retrieve(right_index1);
      result("p,q,i,s") = result("p,q,r,s") * right1("r,i");
    }
    auto right_index2 = formula_string.ket_indices()[1];
    if (right_index2.is_mo()) {
      auto& right2 = orbital_space_registry_->retrieve(right_index2);
      result("p,q,r,i") = result("p,q,r,s") * right2("s,i");
    }

    time1 = mpqc::now(world_, accurate_time_);
    time += mpqc::duration_in_s(time0, time1);
  }

  result.truncate();

  utility::print_par(world_, "Transformed LCAO Integral: ",
                     utility::to_string(formula_string.string()));
  double size = mpqc::detail::array_size(result);
  utility::print_par(world_, " Size: ", size, " GB");
  utility::print_par(world_, " Time: ", time, " s\n");

  return result;
}

template <typename Tile, typename Policy>
Formula LCAOFactory<Tile, Policy>::mo_to_ao(const Formula& formula) {
  std::vector<OrbitalIndex> ao_left_index, ao_right_index;

  int increment = 0;
  auto left_index = formula.bra_indices();
  for (const auto& index : left_index) {
    // find the correspoding ao index
    if (index.is_mo()) {
      auto ao_index = orbital_space_registry_->retrieve(index).ao_index().name();
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
  for (const auto& index : right_index) {
    // find the correspoding ao index
    if (index.is_mo()) {
      auto ao_index = orbital_space_registry_->retrieve(index).ao_index().name();
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

template <typename Tile, typename Policy>
std::tuple<Formula, std::pair<Formula::Position, size_t>>
LCAOFactory<Tile, Policy>::reduce_formula(const Formula& formula) {
  std::vector<float> bra_strength_factors;
  std::vector<float> ket_strength_factors;

  auto compute_strength_factors = [=](
      const std::vector<OrbitalIndex>& indices,
      std::vector<float>& strenth_factors) -> void {
    for (const auto& index : indices) {
      float strength_factor;
      if (index.is_mo()) {
        const auto& orb_space = orbital_space_registry_->retrieve(index);
        const auto rank = orb_space.rank();
        const auto ao_rank = orb_space.ao_rank();
        strength_factor = static_cast<float>(ao_rank) / rank;
      } else
        strength_factor = 0.0;
      strenth_factors.push_back(strength_factor);
    }
  };

  auto bra_indices = formula.bra_indices();
  auto ket_indices = formula.ket_indices();
  TA_USER_ASSERT(!bra_indices.empty() || !ket_indices.empty(),
                 "cannot reduce Formula with empty bra and ket");
  compute_strength_factors(bra_indices, bra_strength_factors);
  compute_strength_factors(ket_indices, ket_strength_factors);

  auto nonzero_float_compare = [](float a, float b) -> bool {
    if (a == 0.0) return false;
    if (b == 0.0) return true;
    return a < b;
  };
  auto min_bra_iter =
      std::min_element(bra_strength_factors.begin(), bra_strength_factors.end(),
                       nonzero_float_compare);
  auto min_ket_iter =
      std::min_element(ket_strength_factors.begin(), ket_strength_factors.end(),
                       nonzero_float_compare);
  Formula::Position pos;  // where the reduced index is located
  if (bra_indices.empty())
    pos = Formula::Position::Ket;
  else if (ket_indices.empty())
    pos = Formula::Position::Bra;
  else {
    pos = *min_bra_iter > *min_ket_iter ? Formula::Position::Bra
                                        : Formula::Position::Ket;
  }

  size_t idx;
  if (pos == Formula::Position::Bra) {
    TA_USER_ASSERT(*min_bra_iter != 0.0,
                   "cannot reduce Formula with AO indices only");
    idx = min_bra_iter - bra_strength_factors.begin();
    auto ao_index =
        orbital_space_registry_->retrieve(bra_indices[idx]).ao_index();
    bra_indices[idx] = ao_index;
  } else {
    TA_USER_ASSERT(*min_ket_iter != 0.0,
                   "cannot reduce Formula with AO indices only");
    idx = min_ket_iter - ket_strength_factors.begin();
    auto ao_index =
        orbital_space_registry_->retrieve(ket_indices[idx]).ao_index();
    ket_indices[idx] = ao_index;
  }

  // reduce the index with maximum strength factor
  auto reduced_formula = formula;
  reduced_formula.set_bra_indices(bra_indices);
  reduced_formula.set_ket_indices(ket_indices);

  return std::make_tuple(reduced_formula, std::make_pair(pos, idx));
}

template <typename Tile, typename Policy>
void LCAOFactory<Tile, Policy>::assert_all_mo(const Formula& formula) {
  auto left = formula.bra_indices();
  for (auto& index : left) {
    TA_ASSERT(index.is_mo());
  }

  auto right = formula.ket_indices();
  for (auto& index : right) {
    TA_ASSERT(index.is_mo());
  }
}
template <typename Tile, typename Policy>
typename LCAOFactory<Tile, Policy>::TArray LCAOFactory<Tile, Policy>::compute(
    const std::wstring& formula_string) {
  Formula formula(formula_string);
  return compute(formula);
}

template <typename Tile, typename Policy>
typename LCAOFactory<Tile, Policy>::TArray LCAOFactory<Tile, Policy>::compute(
    const Formula& formula) {
  auto iter = mo_formula_registry_.find(formula);

  TArray result;

  if (iter != mo_formula_registry_.end()) {
    result = iter->second;
    utility::print_par(world_, "Retrieved LCAO Integral: ",
                       utility::to_string(formula.string()));
    double size = mpqc::detail::array_size(result);
    utility::print_par(world_, " Size: ", size, " GB\n");
    return result;
  } else {
    // find a permutation
    std::vector<Formula> permutes = permutations(formula);
    typename FormulaRegistry<TArray>::iterator find_permute;

    for (auto& permute : permutes) {
      find_permute = mo_formula_registry_.find(permute);
      if (find_permute != mo_formula_registry_.end()) {
        mpqc::time_point time0 = mpqc::now(world_, accurate_time_);

        // permute the array
        result(formula.to_ta_expression()) =
            (find_permute->second)(permute.to_ta_expression());

        mpqc::time_point time1 = mpqc::now(world_, accurate_time_);
        double time = mpqc::duration_in_s(time0, time1);

        utility::print_par(world_, "Permuted LCAO Integral: ",
                           utility::to_string(formula.string()), " From ",
                           utility::to_string(permute.string()));
        double size = mpqc::detail::array_size(result);
        utility::print_par(world_, " Size: ", size, " GB ");
        utility::print_par(world_, " Time: ", time, " s\n");

        // store current array and delete old one
        mo_formula_registry_.insert(formula, result);

        // TODO need to optimize storage and permutation, there is no need to
        // store multiple copy of permutations

        //        mo_formula_registry_.purge_formula(world_,permute);
        return result;
      }
    }

    if (formula.rank() == 2) {
      result = compute2(formula);
      mo_formula_registry_.insert(formula, result);
    } else if (formula.rank() == 3) {
      result = compute3(formula);
      mo_formula_registry_.insert(formula, result);
    } else if (formula.rank() == 4) {
      result = compute4(formula);
      mo_formula_registry_.insert(formula, result);
    }
    madness::print_meminfo(world_.rank(),
                           "LCAOFactory: " + utility::to_string(formula.string()));
    return result;
  }
}

extern template
class LCAOFactory<TA::TensorD,TA::SparsePolicy>;
//extern template
//class LCAOFactory<TA::TensorD,TA::DensePolicy>;

}  // namespace integral
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_LCAO_FACTORY_H_
