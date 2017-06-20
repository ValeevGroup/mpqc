//
// Created by Chong Peng on 1/7/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_LCAO_FACTORY_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_LCAO_FACTORY_H_

#include "mpqc/chemistry/qc/lcao/factory/ao_factory.h"
#include "mpqc/math/linalg/diagonal_array.h"

namespace mpqc {
namespace lcao {

template <typename Tile, typename Policy>
class LCAOFactory;

template <typename Tile, typename Policy>
using LCAOFactoryBase =
    Factory<TA::DistArray<Tile, Policy>,
            TA::DistArray<gaussian::DirectDFTile<Tile, Policy>, Policy>>;

template <typename Tile, typename Policy>
std::shared_ptr<LCAOFactory<Tile, Policy>> construct_lcao_factory(
    const KeyVal& kv) {
  std::shared_ptr<LCAOFactory<Tile, Policy>> lcao_factory;
  if (kv.exists_class("wfn_world:lcao_factory")) {
    lcao_factory =
        kv.class_ptr<LCAOFactory<Tile, Policy>>("wfn_world:lcao_factory");
  } else {
    lcao_factory = std::make_shared<LCAOFactory<Tile, Policy>>(kv);
    std::shared_ptr<DescribedClass> ao_factory_base = lcao_factory;
    KeyVal& kv_nonconst = const_cast<KeyVal&>(kv);
    kv_nonconst.keyval("wfn_world").assign("lcao_factory", ao_factory_base);
  }
  return lcao_factory;
};

template <typename Tile, typename Policy>
LCAOFactory<Tile, Policy>& to_lcao_factory(
    LCAOFactoryBase<Tile, Policy>& factory) {
  return dynamic_cast<LCAOFactory<Tile, Policy>&>(factory);
};

template <typename Tile, typename Policy>
std::shared_ptr<LCAOFactory<Tile, Policy>> to_lcao_factory(
    const std::shared_ptr<LCAOFactoryBase<Tile, Policy>>& factory) {
  auto result = std::dynamic_pointer_cast<LCAOFactory<Tile, Policy>>(factory);
  TA_ASSERT(result != nullptr);
  return result;
};

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
class LCAOFactory : public LCAOFactoryBase<Tile, Policy> {
 public:
  using TArray = TA::DistArray<Tile, Policy>;
  using DirectTArray =
      TA::DistArray<gaussian::DirectDFTile<Tile, Policy>, Policy>;
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
      : LCAOFactoryBase<Tile, Policy>(kv),
        ao_factory_(*gaussian::construct_ao_factory<Tile, Policy>(kv)) {
    std::string prefix = "";
    if (kv.exists("wfn_world") || kv.exists_class("wfn_world")) {
      prefix = "wfn_world:";
    }
    keep_partial_transforms_ =
        kv.value<bool>(prefix + "keep_partial_transform", true);

    auto orbital_space_registry =
        std::make_shared<OrbitalSpaceRegistry<TArray>>();

    this->set_orbital_registry(orbital_space_registry);
    ao_factory_.set_orbital_registry(orbital_space_registry);
    ExEnv::out0() << "\nConstructing LCAOFactory: \n"
                  << indent << "Keep partial transform = "
                  << (keep_partial_transforms_ ? "true" : "false") << "\n"
                  << indent << "Accurate time = "
                  << (this->accurate_time_ ? "true" : "false") << "\n"
                  << indent
                  << "Verbose = " << (this->verbose_ ? "true" : "false")
                  << "\n\n";
  }

  void obsolete() override {
    // obsolete self
    LCAOFactoryBase<Tile, Policy>::obsolete();
    // obsolete AOFactory
    ao_factory_.obsolete();
  }

  /// return reference to AOFactory object
  AOFactoryType& ao_factory() const { return ao_factory_; }

  /// reports the partial tform flag; if true, partially-transformed integrals
  /// are stored
  /// @note at this time only supported for 3-index integrals
  bool keep_partial_transforms() const { return keep_partial_transforms_; }

  /// sets the partial tform flag; if true, partially-transformed integrals are
  /// stored
  /// @note at this time only supported for 3-index integrals
  void keep_partial_transforms(bool flag) { keep_partial_transforms_ = flag; }

  /**
   *  compute integral by Formula
   *  this function will look into registry first
   *  if Formula computed, it will return it from registry
   *  if not, it will compute it
   */
  TArray compute(const Formula&) override;

  DirectTArray compute_direct(const Formula& formula) override;

  using LCAOFactoryBase<Tile, Policy>::compute;
  using LCAOFactoryBase<Tile, Policy>::compute_direct;

  /// purge formula that contain Operator described by string \c str
  /// from mo_registry and ao_registry
  void purge_operator(const std::wstring& str) override {
    LCAOFactoryBase<Tile, Policy>::purge_operator(str);
    ao_factory().purge_operator(str);
  }

  /// purge formulae that contain index described by string \c idx_str
  /// from mo_registry and ao_registry
  void purge_index(const std::wstring& idx_str) override {
    LCAOFactoryBase<Tile, Policy>::purge_index(idx_str);
    ao_factory().purge_index(idx_str);
  }

 private:
  /// compute integrals that has two dimension
  TArray compute2(const Formula& formula_string);

  /// compute integrals that has three dimension
  TArray compute3(const Formula& formula_string);

  /// compute integrals that has four dimension
  TArray compute4(const Formula& formula_string);

  /// "reduce" by converting to AO at most 1 index in the formula,
  /// the index is chosen to minimize strength reduction (thus recursive
  /// reduction ensures maximum strength reduction)
  /// @return {"reduced" formula, {\c pos , \c rank }}, where \c pos and \c rank
  ///         are the location and rank of the reduced index
  std::tuple<Formula, std::pair<Formula::Position, size_t>> reduce_formula(
      const Formula& formula);

 private:
  AOFactoryType& ao_factory_;

  bool keep_partial_transforms_;  //!< if true, keep partially-transformed ints
                                  //!(false by default)
};

template <typename Tile, typename Policy>
typename LCAOFactory<Tile, Policy>::TArray LCAOFactory<Tile, Policy>::compute2(
    const Formula& formula_string) {
  double time = 0.0;
  mpqc::time_point time0;
  mpqc::time_point time1;
  auto& world = this->world();

  TArray result;
  // Identity matrix
  if (formula_string.oper().type() == Operator::Type::Identity) {
    time0 = mpqc::now(world, this->accurate_time_);

    auto left_index1 = formula_string.bra_indices()[0];
    auto right_index1 = formula_string.ket_indices()[0];
    auto& left1 = this->orbital_registry().retrieve(left_index1);
    auto& right1 = this->orbital_registry().retrieve(right_index1);

    auto tr1 = left1.trange();
    auto tr2 = right1.trange();

    // create diagonal array
    result = array_ops::create_diagonal_array_from_eigen<Tile, Policy>(
        world, tr1, tr2, 1.0);
    result.truncate();

    time1 = mpqc::now(world, this->accurate_time_);
    time += mpqc::duration_in_s(time0, time1);

    if (this->verbose_) {
      ExEnv::out0() << indent;
      ExEnv::out0() << "Computed Identity: "
                    << utility::to_string(formula_string.string());
      double size = mpqc::detail::array_size(result);
      ExEnv::out0() << " Size: " << size << " GB"
                    << " Time: " << time << " s\n";
    }
    return result;
  }

  // get AO
  auto ao_formula =
      detail::lcao_to_ao(formula_string, this->orbital_registry());
  auto ao_factory = ao_factory_.compute(ao_formula);

  time0 = mpqc::now(world, this->accurate_time_);
  // convert to MO
  result = ao_factory;
  // get coefficient
  auto left_index1 = formula_string.bra_indices()[0];
  if (left_index1.is_lcao()) {
    auto& left1 = this->orbital_registry().retrieve(left_index1);
    result("i,r") = result("p,r") * left1("p,i");
  }
  auto right_index1 = formula_string.ket_indices()[0];
  if (right_index1.is_lcao()) {
    auto& right1 = this->orbital_registry().retrieve(right_index1);
    result("p,k") = result("p,r") * right1("r,k");
  }

  result.truncate();
  time1 = mpqc::now(world, this->accurate_time_);
  time += mpqc::duration_in_s(time0, time1);

  if (this->verbose_) {
    ExEnv::out0() << "Transformed LCAO Integral: "
                  << utility::to_string(formula_string.string());
    double size = mpqc::detail::array_size(result);
    ExEnv::out0() << " Size: " << size << " GB"
                  << " Time: " << time << " s\n";
  }
  return result;
}

template <typename Tile, typename Policy>
typename LCAOFactory<Tile, Policy>::TArray LCAOFactory<Tile, Policy>::compute3(
    const Formula& formula_string) {
  double time = 0.0;
  mpqc::time_point time0;
  mpqc::time_point time1;
  auto& world = this->world();

  TArray result;

  if (formula_string.has_option(Formula::Option::Inverse) ||
      formula_string.has_option(Formula::Option::InverseSquareRoot)) {
    auto three_center = formula_string;
    three_center.clear_option();
    TArray three_center_array;

    auto iter = this->registry_.find(three_center);
    // if three center array already in registry, use it
    if( iter != this->registry_.end() ){
      three_center_array = iter->second;
    }
    else{
      // do not store temporary three center array in registry
      three_center_array = this->compute3(three_center);
    }

    auto two_center = formula_string;
    two_center.set_ket_indices(formula_string.bra_indices());
    auto two_center_array = this->ao_factory_.compute(two_center);

    time0 = mpqc::now(world, this->accurate_time_);
    result("K,i,j") = two_center_array("K,Q") * three_center_array("Q,i,j");

  } else {
    if (not keep_partial_transforms()) {  // compute from AO ints
      // get AO
      auto ao_formula =
          detail::lcao_to_ao(formula_string, this->orbital_registry());
      auto ao_integral = ao_factory_.compute_direct(ao_formula);

      time0 = mpqc::now(world, this->accurate_time_);

      // transform to MO, only convert the right side
      auto right_index1 = formula_string.ket_indices()[0];
      if (right_index1.is_lcao()) {
        auto& right1 = this->orbital_registry().retrieve(right_index1);
        result("K,i,q") = ao_integral("K,p,q") * right1("p,i");
      }
      auto right_index2 = formula_string.ket_indices()[1];
      if (right_index2.is_lcao()) {
        auto& right2 = this->orbital_registry().retrieve(right_index2);
        result("K,p,j") = result("K,p,q") * right2("q,j");
      }
    } else {  // tform to optimally reduce strength, store partial transform
              // results

      // compute reduced formula
      Formula reduced_formula;  // reduced formula
      std::pair<Formula::Position, size_t> reduced_index_coord;
      std::tie(reduced_formula, reduced_index_coord) =
          reduce_formula(formula_string);

      const auto reduced_index_position = reduced_index_coord.first;
      const auto reduced_index_rank = reduced_index_coord.second;
      const auto reduced_index_absrank =
          reduced_index_rank + (reduced_index_position == Formula::Position::Bra
                                    ? 0
                                    : reduced_formula.bra_indices().size());

      // extract orbital coefficients for transformation
      auto reduced_index =
          (reduced_index_position == Formula::Position::Bra)
              ? formula_string.bra_indices()[reduced_index_rank]
              : formula_string.ket_indices()[reduced_index_rank];
      auto& reduced_index_space =
          this->orbital_registry().retrieve(reduced_index);
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

      // compute direct three center AO
      if (reduced_formula.is_ao() &&
          !ao_factory_.registry().have(reduced_formula)) {
        auto reduced_integral = ao_factory_.compute_direct(reduced_formula);
        time0 = mpqc::now(world, this->accurate_time_);
        result(result_key) =
            reduced_integral(reduced_key) * reduced_index_coeff(coeff_key);
      } else {
        auto reduced_integral = reduced_formula.is_ao()
                                    ? ao_factory_.compute(reduced_formula)
                                    : this->compute(reduced_formula);
        time0 = mpqc::now(world, this->accurate_time_);
        result(result_key) =
            reduced_integral(reduced_key) * reduced_index_coeff(coeff_key);
      }
    }
  }

  result.truncate();

  time1 = mpqc::now(world, this->accurate_time_);
  time += mpqc::duration_in_s(time0, time1);

  if (this->verbose_) {
    ExEnv::out0() << indent;
    ExEnv::out0() << "Transformed LCAO Integral: "
                  << utility::to_string(formula_string.string());
    double size = mpqc::detail::array_size(result);
    ExEnv::out0() << " Size: " << size << " GB"
                  << " Time: " << time << " s\n";
  }

  return result;
};

template <typename Tile, typename Policy>
typename LCAOFactory<Tile, Policy>::TArray LCAOFactory<Tile, Policy>::compute4(
    const Formula& formula_string) {
  double time = 0.0;
  mpqc::time_point time0;
  mpqc::time_point time1;
  auto& world = this->world();
  TArray result;
  if (formula_string.has_option(Formula::Option::DensityFitting)) {
    // get df formula
    auto df_formulas = gaussian::detail::get_df_formula(formula_string);

    // compute integral
    TArray left = compute(df_formulas[0]);

    TArray right = compute(df_formulas[1]);

    time0 = mpqc::now(world, this->accurate_time_);

    auto notation = formula_string.notation();
    std::string result_str =
        notation == Formula::Notation::Chemical ? "i,j,k,l" : "i,k,j,l";

    if (formula_string.oper().type() == Operator::Type::cGTG ||
        formula_string.oper().type() == Operator::Type::cGTGCoulomb) {
      result(result_str) = left("q,i,j") * (-right("q,k,l"));
    } else {
      result(result_str) = left("q,i,j") * right("q,k,l");
    }

    time1 = mpqc::now(world, this->accurate_time_);
    time += mpqc::duration_in_s(time0, time1);

  } else {
    // get AO
    auto ao_formula =
        detail::lcao_to_ao(formula_string, this->orbital_registry());

    // convert to MO
    time0 = mpqc::now(world, this->accurate_time_);

    // get coefficient
    auto left_index1 = formula_string.bra_indices()[0];

    // if first index is occupied, store half transformed integral
    TA_ASSERT(left_index1.is_lcao());
    if (left_index1 == OrbitalIndex(L"i") ||
        left_index1 == OrbitalIndex(L"m")) {
      // make ( i rho | sigma mu )
      auto half_transformed_formula = ao_formula;
      auto bra_ind = half_transformed_formula.bra_indices();
      bra_ind[0] = left_index1;
      half_transformed_formula.set_bra_indices(bra_ind);
      if (half_transformed_formula.notation() == Formula::Notation::Physical) {
        half_transformed_formula.set_notation(Formula::Notation::Chemical);
        ao_formula.set_notation(Formula::Notation::Chemical);
      }

      auto iter = this->registry_.find(half_transformed_formula);
      if (iter != this->registry_.end()) {
        result = this->compute(half_transformed_formula);
      } else {
        auto direct_ao = ao_factory_.compute_direct(ao_formula);
        auto& left1 = this->orbital_registry().retrieve(left_index1);
        result("i, rho, mu, nu") =
            left1("sigma, i") * direct_ao("sigma, rho, mu, nu");

        this->registry_.insert(half_transformed_formula, result);
      }
    }
    // if first index not occupied, intermediate too large, not store
    // intermediate
    else {
      auto half_transformed_formula = ao_formula;
      auto bra_ind = half_transformed_formula.bra_indices();
      bra_ind[0] = left_index1;
      half_transformed_formula.set_bra_indices(bra_ind);
      if (half_transformed_formula.notation() == Formula::Notation::Physical) {
        half_transformed_formula.set_notation(Formula::Notation::Chemical);
        ao_formula.set_notation(Formula::Notation::Chemical);
      }
      if (ao_formula.notation() == Formula::Notation::Physical) {
        ao_formula.set_notation(Formula::Notation::Chemical);
      }
      auto direct_ao = ao_factory_.compute_direct(ao_formula);
      auto& left1 = this->orbital_registry().retrieve(left_index1);
      result("a, rho, mu, nu") =
          left1("sigma, a") * direct_ao("sigma, rho, mu, nu");
      ExEnv::out0() << indent;
      ExEnv::out0() << "Waring! Transformation creates large intermediate:  ";
      ExEnv::out0() << utility::to_string(half_transformed_formula.string());
      double size = mpqc::detail::array_size(result);
      ExEnv::out0() << " Size: " << size << " GB\n";
    }

    auto left_index2 = formula_string.bra_indices()[1];
    if (left_index2.is_lcao()) {
      auto& left2 = this->orbital_registry().retrieve(left_index2);
      if (formula_string.notation() == Formula::Notation::Physical) {
        result("p,i,q,s") = result("p,q,r,s") * left2("r,i");
      } else {
        result("p,i,r,s") = result("p,q,r,s") * left2("q,i");
      }
    }

    auto right_index1 = formula_string.ket_indices()[0];
    if (right_index1.is_lcao()) {
      auto& right1 = this->orbital_registry().retrieve(right_index1);
      result("p,q,i,s") = result("p,q,r,s") * right1("r,i");
    }
    auto right_index2 = formula_string.ket_indices()[1];
    if (right_index2.is_lcao()) {
      auto& right2 = this->orbital_registry().retrieve(right_index2);
      result("p,q,r,i") = result("p,q,r,s") * right2("s,i");
    }

    time1 = mpqc::now(world, this->accurate_time_);
    time += mpqc::duration_in_s(time0, time1);
  }

  result.truncate();

  if (this->verbose_) {
    ExEnv::out0() << indent;
    ExEnv::out0() << "Transformed LCAO Integral: "
                  << utility::to_string(formula_string.string());
    double size = mpqc::detail::array_size(result);
    ExEnv::out0() << " Size: " << size << " GB"
                  << " Time: " << time << " s\n";
  }

  return result;
}

template <typename Tile, typename Policy>
std::tuple<Formula, std::pair<Formula::Position, size_t>>
LCAOFactory<Tile, Policy>::reduce_formula(const Formula& formula) {
  std::vector<float> bra_strength_factors;
  std::vector<float> ket_strength_factors;

  auto compute_strength_factors = [=](const std::vector<OrbitalIndex>& indices,
                                      std::vector<float>& strenth_factors) {
    for (const auto& index : indices) {
      float strength_factor;
      if (index.is_lcao()) {
        const auto& orb_space = this->orbital_registry().retrieve(index);
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

  auto nonzero_float_compare = [](float a, float b) {
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
        this->orbital_registry().retrieve(bra_indices[idx]).ao_index();
    bra_indices[idx] = ao_index;
  } else {
    TA_USER_ASSERT(*min_ket_iter != 0.0,
                   "cannot reduce Formula with AO indices only");
    idx = min_ket_iter - ket_strength_factors.begin();
    auto ao_index =
        this->orbital_registry().retrieve(ket_indices[idx]).ao_index();
    ket_indices[idx] = ao_index;
  }

  // reduce the index with maximum strength factor
  auto reduced_formula = formula;
  reduced_formula.set_bra_indices(bra_indices);
  reduced_formula.set_ket_indices(ket_indices);

  return std::make_tuple(reduced_formula, std::make_pair(pos, idx));
}

template <typename Tile, typename Policy>
typename LCAOFactory<Tile, Policy>::TArray LCAOFactory<Tile, Policy>::compute(
    const Formula& formula) {
  TA_USER_ASSERT(!lcao::detail::if_all_ao(formula),
                 "LCAOFactory cannot accept all AO index!\n");
  ExEnv::out0() << incindent;

  auto& world = this->world();

  auto iter = this->registry_.find(formula);

  TArray result;

  if (iter != this->registry_.end()) {
    result = iter->second;

    if (this->verbose_) {
      ExEnv::out0() << indent;
      ExEnv::out0() << "Retrieved LCAO Integral: "
                    << utility::to_string(formula.string());
      double size = mpqc::detail::array_size(result);
      ExEnv::out0() << " Size: " << size << " GB\n";
    }
  } else {
    // find a permutation
    std::vector<Formula> permutes = permutations(formula);
    typename FormulaRegistry<TArray>::iterator find_permute;

    for (auto& permute : permutes) {
      find_permute = this->registry_.find(permute);
      if (find_permute != this->registry_.end()) {
        mpqc::time_point time0 = mpqc::now(world, this->accurate_time_);

        // permute the array
        result(formula.to_ta_expression()) =
            (find_permute->second)(permute.to_ta_expression());

        mpqc::time_point time1 = mpqc::now(world, this->accurate_time_);
        double time = mpqc::duration_in_s(time0, time1);

        if (this->verbose_) {
          ExEnv::out0() << indent;
          ExEnv::out0() << "Permuted LCAO Integral: "
                        << utility::to_string(formula.string()) << " From "
                        << utility::to_string(permute.string());
          double size = mpqc::detail::array_size(result);
          ExEnv::out0() << " Size: " << size << " GB "
                        << " Time: " << time << " s\n";
        }

        // store current array and delete old one
        this->registry_.insert(formula, result);
        this->registry_.purge_formula(permute);
      }
    }

    // if not find formula

    if (!result.is_initialized()) {
      if (formula.rank() == 2) {
        result = compute2(formula);
        this->registry_.insert(formula, result);
      } else if (formula.rank() == 3) {
        result = compute3(formula);
        this->registry_.insert(formula, result);
      } else if (formula.rank() == 4) {
        result = compute4(formula);
        this->registry_.insert(formula, result);
      }
    }

    madness::print_meminfo(
        world.rank(), "LCAOFactory: " + utility::to_string(formula.string()));
  }
  ExEnv::out0() << decindent;
  return result;
}

template <typename Tile, typename Policy>
typename LCAOFactory<Tile, Policy>::DirectTArray
LCAOFactory<Tile, Policy>::compute_direct(const Formula& formula) {
  TA_ASSERT(formula.rank() == 4);
  TA_ASSERT(formula.has_option(Formula::Option::DensityFitting));

  ExEnv::out0() << incindent;
  double time = 0.0;
  mpqc::time_point time0;
  mpqc::time_point time1;
  auto& world = this->world();
  DirectTArray result;

  auto iter = this->direct_registry_.find(formula);

  if (iter != this->direct_registry_.end()) {
    result = iter->second;

    if (this->verbose_) {
      ExEnv::out0() << indent;
      ExEnv::out0() << "Retrieved LCAO Direct Integral From Density-Fitting: "
                    << utility::to_string(formula.string());
      double size = mpqc::detail::array_size(result);
      ExEnv::out0() << " Size: " << size << " GB\n";
    }
  } else {
    // get three center integral
    auto df_formulas = gaussian::detail::get_df_formula(formula);
    TArray left = compute(df_formulas[0]);
    TArray right = compute(df_formulas[1]);

    time0 = mpqc::now(world, this->accurate_time_);

    result = gaussian::df_direct_integrals(left, right, formula.notation());

    time1 = mpqc::now(world, this->accurate_time_);
    time = mpqc::duration_in_s(time0, time1);
    this->direct_registry_.insert(formula, result);
  }

  if (this->verbose_) {
    ExEnv::out0() << indent;
    ExEnv::out0() << "Computed LCAO Direct Integral From Density-Fitting: "
                  << utility::to_string(formula.string());
    double size = mpqc::detail::array_size(result);
    ExEnv::out0() << " Size: " << size << " GB"
                  << " Time: " << time << " s\n";
    ExEnv::out0() << decindent;
  }
  return result;
};

#if TA_DEFAULT_POLICY == 0
extern template class LCAOFactory<TA::TensorD, TA::DensePolicy>;
#elif TA_DEFAULT_POLICY == 1
extern template class LCAOFactory<TA::TensorD, TA::SparsePolicy>;
#endif

}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_LCAO_FACTORY_H_
