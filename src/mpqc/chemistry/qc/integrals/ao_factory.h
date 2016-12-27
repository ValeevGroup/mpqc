//
// Created by Chong Peng on 10/14/15.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_AO_FACTORY_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_AO_FACTORY_H_

#include "ao_factory_base.h"
#include "mpqc/chemistry/qc/expression/permutation.h"
#include "mpqc/chemistry/qc/integrals/f12_utility.h"
#include "mpqc/chemistry/qc/integrals/integrals.h"
#include "mpqc/math/external/eigen/eigen.h"
#include "mpqc/math/linalg/sqrt_inv.h"
#include "mpqc/util/external/madworld/parallel_break_point.h"
#include "mpqc/util/external/madworld/parallel_print.h"
#include "mpqc/util/misc/time.h"
#include <madness/world/worldmem.h>

namespace mpqc {
namespace lcao {
namespace gaussian {

template <typename Tile, typename Policy>
class AOFactory;

template <typename Tile, typename Policy>
std::shared_ptr<AOFactory<Tile, Policy>> construct_ao_factory(
    const KeyVal& kv) {
  std::shared_ptr<AOFactory<Tile, Policy>> ao_factory;
  if (kv.exists_class("wfn_world:ao_factory")) {
    ao_factory = kv.class_ptr<AOFactory<Tile, Policy>>("wfn_world:ao_factory");
  } else {
    ao_factory = std::make_shared<AOFactory<Tile, Policy>>(kv);
    std::shared_ptr<DescribedClass> ao_factory_base = ao_factory;
    KeyVal& kv_nonconst = const_cast<KeyVal&>(kv);
    kv_nonconst.keyval("wfn_world").assign("ao_factory", ao_factory_base);
  }
  return ao_factory;
};

// TODO better inverse of two center
// TODO direct integral
// TODO Screener for different type of integral
/**
 * \brief Atomic Integral Class
 *
 *  This class computes atomic integral using  Formula
 *
 *  compute(formula) return TArray object
 *  (formula) return TArray expression
 *
 */

template <typename Tile, typename Policy = TA::SparsePolicy>
class AOFactory : public AOFactoryBase, public DescribedClass {
 public:
  using TArray = TA::DistArray<Tile, Policy>;

  /// Op is a function pointer that convert TA::Tensor to Tile
  //  using Op = Tile (*) (TA::TensorD&&);
  using Op = std::function<Tile(TA::TensorD&&)>;

  AOFactory() = default;

  /**
   * \brief  KeyVal constructor
   *
   * It takes all the keys to construct AOFactoryBase and also the following
   *
   *  | KeyWord | Type | Default| Description |
   *  |---------|------|--------|-------------|
   *  |accurate_time|bool|false|if do fence at timing|
   *  |iterative_inv_sqrt|bool|false| use iterative inverse square root |
   *
   *
   **/

  AOFactory(const KeyVal& kv)
      : AOFactoryBase(kv), ao_formula_registry_(), orbital_space_registry_() {
    std::string prefix = "";
    if (kv.exists("wfn_world") || kv.exists_class("wfn_world")) {
      prefix = "wfn_world:";
    }

    accurate_time_ = kv.value<bool>(prefix + "accurate_time", false);
    iterative_inv_sqrt_ = kv.value<bool>(prefix + "iterative_inv_sqrt", false);

    set_oper(Tile());
  }

  AOFactory(AOFactory&&) = default;
  AOFactory& operator=(AOFactory&&) = default;

  /// set oper based on Tile type
  template <typename T = Tile>
  void set_oper(typename std::enable_if<std::is_same<T, TA::TensorD>::value,
                                        T>::type&& t) {
    op_ = TA::Noop<TA::TensorD, true>();
  }

  virtual ~AOFactory() noexcept = default;

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
    auto array = compute(formula);
    auto& result = ao_formula_registry_.retrieve(formula);
    return result(formula.to_ta_expression());
  };

  /// return ao formula registry
  const FormulaRegistry<TArray>& registry() const {
    return ao_formula_registry_;
  }

  /// return ao formula registry
  FormulaRegistry<TArray>& registry() { return ao_formula_registry_; }

  /// set orbital space registry
  void set_orbital_space_registry(
      const std::shared_ptr<OrbitalSpaceRegistry<TArray>> regitsry) {
    orbital_space_registry_ = regitsry;
  }

 protected:
  /// compute integrals that has two dimension
  TArray compute2(const Formula& formula_string);

  /// compute integrals that has three dimension
  TArray compute3(const Formula& formula_string);

  /// compute integrals that has four dimension
  TArray compute4(const Formula& formula_string);

 private:
  /// compute sparse array
  template <typename U = Policy>
  TA::DistArray<
      Tile, typename std::enable_if<std::is_same<U, TA::SparsePolicy>::value,
                                    TA::SparsePolicy>::type>
  compute_integrals(
      madness::World& world, ShrPool<libint2::Engine>& engine,
      BasisVector const& bases,
      std::shared_ptr<Screener> p_screen =
          std::make_shared<Screener>(Screener{})) {
    auto result =
        sparse_integrals(world, engine, bases, p_screen, op_);
    return result;
  }

  /// compute dense array
  template <typename U = Policy>
  TA::DistArray<Tile,
                typename std::enable_if<std::is_same<U, TA::DensePolicy>::value,
                                        TA::DensePolicy>::type>
  compute_integrals(
      madness::World& world, ShrPool<libint2::Engine>& engine,
      BasisVector const& bases,
      std::shared_ptr<Screener> p_screen =
          std::make_shared<Screener>(Screener{})) {
    auto result =
        dense_integrals(world, engine, bases, p_screen, op_);
    return result;
  }

 private:
  FormulaRegistry<TArray> ao_formula_registry_;
  std::shared_ptr<OrbitalSpaceRegistry<TArray>> orbital_space_registry_;
  Op op_;
  bool accurate_time_;
  bool iterative_inv_sqrt_;
};

template <typename Tile, typename Policy>
typename AOFactory<Tile, Policy>::TArray AOFactory<Tile, Policy>::compute(
    const std::wstring& formula_string) {
  auto formula = Formula(formula_string);
  return compute(formula);
}

template <typename Tile, typename Policy>
typename AOFactory<Tile, Policy>::TArray AOFactory<Tile, Policy>::compute(
    const Formula& formula) {
  TArray result;

  // find if in registry
  auto iter = ao_formula_registry_.find(formula);

  if (iter != ao_formula_registry_.end()) {
    result = *(iter->second);
    utility::print_par(world_, "Retrieved AO Integral: ",
                       utility::to_string(formula.string()));
    double size = mpqc::detail::array_size(result);
    utility::print_par(world_, " Size: ", size, " GB\n");
    return result;
  } else {
    // find a permutation, currently, it won't store permutation in registry

    std::vector<Formula> permutes = permutations(formula);
    typename FormulaRegistry<TArray>::iterator find_permute;

    for (auto& permute : permutes) {
      find_permute = ao_formula_registry_.find(permute);
      if (find_permute != ao_formula_registry_.end()) {
        mpqc::time_point time0 = mpqc::now(world_, accurate_time_);

        // permute the array
        result(formula.to_ta_expression()) =
            (*(find_permute->second))(permute.to_ta_expression());

        mpqc::time_point time1 = mpqc::now(world_, accurate_time_);
        double time = mpqc::duration_in_s(time0, time1);

        utility::print_par(world_, "Permuted AO Integral: ",
                           utility::to_string(formula.string()), " From ",
                           utility::to_string(permute.string()));
        double size = mpqc::detail::array_size(result);
        utility::print_par(world_, " Size: ", size, " GB ");
        utility::print_par(world_, " Time: ", time, " s\n");

        // store current array and delete old one
        ao_formula_registry_.insert(formula, result);
        //        ao_formula_registry_.purge_formula(world_,permute);
        return result;
      }
    }

    // compute formula
    if (formula.rank() == 2) {
      result = compute2(formula);
      ao_formula_registry_.insert(formula, result);
    } else if (formula.rank() == 3) {
      result = compute3(formula);
      ao_formula_registry_.insert(formula, result);
    } else if (formula.rank() == 4) {
      result = compute4(formula);
      ao_formula_registry_.insert(formula, result);
    }

    madness::print_meminfo(
        world_.rank(), "AOFactory: " + utility::to_string(formula.string()));

    return ao_formula_registry_.retrieve(formula);
  }
  // wait all process to obtain and insert result
  //        world_.gop.fence();
}

template <typename Tile, typename Policy>
typename AOFactory<Tile, Policy>::TArray AOFactory<Tile, Policy>::compute2(
    const Formula& formula) {
  BasisVector bs_array;
  double time = 0.0;
  mpqc::time_point time0;
  mpqc::time_point time1;
  TArray result;

  // get the inverse square root instead
  if (iterative_inv_sqrt_ &&
      formula.oper().has_option(Operator::Option::Inverse)) {
    auto inv_sqrt_formula = formula;
    inv_sqrt_formula.set_operator_option({Operator::Option::InverseSquareRoot});

    result = this->compute(inv_sqrt_formula);

    time0 = mpqc::now(world_, accurate_time_);
    result("p,q") = result("p,r") * result("r,q");

    if (formula.oper().type() == Operator::Type::cGTG ||
        formula.oper().type() == Operator::Type::cGTGCoulomb) {
      result("i,j") = -result("i,j");
    }
    utility::print_par(world_, "Computed Inverse of Integral: ",
                       utility::to_string(formula.string()));
    double size = mpqc::detail::array_size(result);

    time1 = mpqc::now(world_, accurate_time_);
    time += mpqc::duration_in_s(time0, time1);
    utility::print_par(world_, " Size: ", size, " GB");
    utility::print_par(world_, " Time: ", time, " s\n");
  }
  // continue with normal step
  else {
    // use one body engine
    if (formula.oper().is_onebody()) {
      // H = V + T
      if (formula.oper().type() == Operator::Type::Core) {
        auto v_formula = formula;
        v_formula.set_operator_type(Operator::Type::Nuclear);

        auto t_formula = formula;
        t_formula.set_operator_type(Operator::Type::Kinetic);

        auto v = this->compute(v_formula);
        auto t = this->compute(t_formula);

        time0 = mpqc::now(world_, accurate_time_);

        result("i,j") = v("i,j") + t("i,j");

        time1 = mpqc::now(world_, accurate_time_);
        time += mpqc::duration_in_s(time0, time1);
      }
      // one body integral S, V, T...
      else {
        time0 = mpqc::now(world_, accurate_time_);

        std::shared_ptr<utility::TSPool<libint2::Engine>> engine_pool;
        parse_one_body(formula, engine_pool, bs_array);
        result = compute_integrals(this->world_, engine_pool, bs_array);

        time1 = mpqc::now(world_, accurate_time_);
        time += mpqc::duration_in_s(time0, time1);
      }

      utility::print_par(world_, "Computed One Body Integral: ",
                         utility::to_string(formula.string()));
      double size = mpqc::detail::array_size(result);
      utility::print_par(world_, " Size: ", size, " GB");
      utility::print_par(world_, " Time: ", time, " s\n");
    }
    // use two body engine
    else if (formula.oper().is_twobody()) {
      time0 = mpqc::now(world_, accurate_time_);

      // compute integral
      std::shared_ptr<utility::TSPool<libint2::Engine>> engine_pool;
      parse_two_body_two_center(formula, engine_pool, bs_array);
      result = compute_integrals(this->world_, engine_pool, bs_array);


      time1 = mpqc::now(world_, accurate_time_);
      time += mpqc::duration_in_s(time0, time1);

      utility::print_par(world_, "Computed Twobody Two Center Integral: ",
                         utility::to_string(formula.string()));
      double size = mpqc::detail::array_size(result);
      utility::print_par(world_, " Size: ", size, " GB");
      utility::print_par(world_, " Time: ", time, " s\n");
    }
    // compute JK, requires orbital space registry
    else if (formula.oper().is_jk()) {
      // density fitting case

      // find the density
      auto space_index = get_jk_orbital_space(formula.oper());
      auto& space = orbital_space_registry_->retrieve(space_index);

      auto obs = space.ao_index().name();
      if (formula.oper().has_option(Operator::Option::DensityFitting)) {
        auto three_center_formula = get_jk_df_formula(formula, obs);

        auto left = compute(three_center_formula[0]);
        auto center = compute(three_center_formula[1]);
        auto right = compute(three_center_formula[2]);

        time0 = mpqc::now(world_, accurate_time_);

        // J case
        if (formula.oper().type() == Operator::Type::J) {
          result("i,j") = center("K,Q") * right("Q,k,l") *
                          (space("k,a") * space("l,a")) * left("K,i,j");
        }
        // K case
        else {
          result("i,j") = (left("K,i,k") * space("k,a")) * center("K,Q") *
                          (right("Q,j,l") * space("l,a"));
        }

        time1 = mpqc::now(world_, accurate_time_);
        time += mpqc::duration_in_s(time0, time1);
      }
      // four center case
      else {
        // find the density
        auto space_index = get_jk_orbital_space(formula.oper());
        auto& space = orbital_space_registry_->retrieve(space_index);
        auto obs = space.ao_index().name();
        // convert to ao formula
        auto four_center_formula = get_jk_formula(formula, obs);
        auto four_center = this->compute(four_center_formula);

        time0 = mpqc::now(world_, accurate_time_);

        if (formula.notation() == Formula::Notation::Chemical) {
          if (formula.oper().type() == Operator::Type::J) {
            result("rho,sigma") = four_center("rho,sigma,mu,nu") *
                                  (space("mu,i") * space("nu,i"));
          } else {
            result("rho,sigma") = four_center("rho,mu,sigma,nu") *
                                  (space("mu,i") * space("nu,i"));
          }
        } else {
          if (formula.oper().type() == Operator::Type::K) {
            result("rho,sigma") = four_center("rho,sigma,mu,nu") *
                                  (space("mu,i") * space("nu,i"));
          } else {
            result("rho,sigma") = four_center("rho,mu,sigma,nu") *
                                  (space("mu,i") * space("nu,i"));
          }
        }

        time1 = mpqc::now(world_, accurate_time_);
        time += mpqc::duration_in_s(time0, time1);
      }
      utility::print_par(world_, "Computed Coulumb/Exchange Integral: ",
                         utility::to_string(formula.string()));
      double size = mpqc::detail::array_size(result);
      utility::print_par(world_, " Size: ", size, " GB");
      utility::print_par(world_, " Time: ", time, " s\n");
    }
    // hJ = H + J
    else if (formula.oper().type() == Operator::Type::hJ) {
      auto h_formula = formula;
      h_formula.set_operator(Operator(L"H"));

      auto j_formula = formula;
      j_formula.set_operator_type(Operator::Type::J);

      auto h = this->compute(h_formula);
      auto j = this->compute(j_formula);

      time0 = mpqc::now(world_, accurate_time_);

      result("i,j") = h("i,j") + 2 * j("i,j");

      time1 = mpqc::now(world_, accurate_time_);
      time += mpqc::duration_in_s(time0, time1);

      utility::print_par(world_, "Computed Coulumb/Exchange Integral: ",
                         utility::to_string(formula.string()));
      double size = mpqc::detail::array_size(result);
      utility::print_par(world_, " Size: ", size, " GB");
      utility::print_par(world_, " Time: ", time, " s\n");
    }
    // compute Fock, requires orbital space registry
    else if (formula.oper().is_fock()) {
      auto formulas = get_fock_formula(formula);

      auto h = compute(formulas[0]);
      auto j = compute(formulas[1]);
      auto k = compute(formulas[2]);

      time0 = mpqc::now(world_, accurate_time_);
      // if closed shell
      if (formula.oper().type() == Operator::Type::Fock) {
        result("rho,sigma") =
            h("rho,sigma") + 2 * j("rho,sigma") - k("rho,sigma");
      }
      // else if spin orbital
      else {
        result("rho,sigma") = h("rho,sigma") + j("rho,sigma") - k("rho,sigma");
      }

      time1 = mpqc::now(world_, accurate_time_);

      time += mpqc::duration_in_s(time0, time1);

      utility::print_par(world_, "Computed Fock Integral: ",
                         utility::to_string(formula.string()));
      double size = mpqc::detail::array_size(result);
      utility::print_par(world_, " Size: ", size, " GB");
      utility::print_par(world_, " Time: ", time, " s\n");
    }
  }

  // check for other options

  // compute inverse square root first in this case
  if (!iterative_inv_sqrt_ &&
      formula.oper().has_option(Operator::Option::Inverse)) {

    time0 = mpqc::now(world_, accurate_time_);

    if (formula.oper().type() == Operator::Type::cGTG ||
        formula.oper().type() == Operator::Type::cGTGCoulomb) {
      result("i,j") = -result("i,j");
    }

    RowMatrixXd result_eig = array_ops::array_to_eigen(result);

    // compute cholesky decomposition
    auto llt_solver = Eigen::LLT<RowMatrixXd>(result_eig);

    // check success
    Eigen::ComputationInfo info = llt_solver.info();
    if (info == Eigen::ComputationInfo::Success) {
      RowMatrixXd L = RowMatrixXd(llt_solver.matrixL());
      RowMatrixXd L_inv_eig = L.inverse();
      result_eig = L_inv_eig.transpose() * L_inv_eig;
    } else if (info == Eigen::ComputationInfo::NumericalIssue) {
      utility::print_par(world_,
                         "!!!\nWarning!! NumericalIssue in Cholesky "
                             "Decomposition\n!!!\n");
    } else if (info == Eigen::ComputationInfo::NoConvergence) {
      utility::print_par(world_,
                         "!!!\nWarning!! NoConvergence in Cholesky "
                             "Decomposition\n!!!\n");
    }

    if (info != Eigen::ComputationInfo::Success) {
      utility::print_par(world_, "Using Eigen LU Decomposition Inverse!\n");

      Eigen::FullPivLU<RowMatrixXd> lu(result_eig);

      TA_ASSERT(lu.isInvertible());

      result_eig = lu.inverse();
    }

    auto tr_result = result.trange().data()[0];
    result = array_ops::eigen_to_array<Tile,Policy>(
        result.world(), result_eig, tr_result, tr_result);

    if (formula.oper().type() == Operator::Type::cGTG ||
        formula.oper().type() == Operator::Type::cGTGCoulomb) {
      result("i,j") = -result("i,j");
    }
    time1 = mpqc::now(world_, accurate_time_);
    auto inv_time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world_, "Inverse Time: ", inv_time, " s\n");
  }

  // inverse square root of integral
  if (formula.oper().has_option(Operator::Option::InverseSquareRoot)) {
    time0 = mpqc::now(world_, accurate_time_);
    if (formula.oper().type() == Operator::Type::cGTG ||
        formula.oper().type() == Operator::Type::cGTGCoulomb) {
      result("i,j") = -result("i,j");
    }

    if (iterative_inv_sqrt_) {
      TArray tmp = array_ops::inverse_sqrt(result);
      result("i,j") = tmp("i,j");
    } else {
      auto result_eig = array_ops::array_to_eigen(result);

      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(result_eig);
      RowMatrixXd inv_eig = es.operatorInverseSqrt();

      auto tr_result = result.trange().data()[0];
      result = array_ops::eigen_to_array<Tile,Policy>(
          result.world(), inv_eig, tr_result, tr_result);
    }

    if (formula.oper().type() == Operator::Type::cGTG ||
        formula.oper().type() == Operator::Type::cGTGCoulomb) {
      result("i,j") = -result("i,j");
    }
    time1 = mpqc::now(world_, accurate_time_);
    auto inv_sqrt_time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world_, "Inverse Square Root Time: ", inv_sqrt_time, " s\n");
  }

  return result;
}

template <typename Tile, typename Policy>
typename AOFactory<Tile, Policy>::TArray AOFactory<Tile, Policy>::compute3(
    const Formula& formula) {
  double time = 0.0;
  mpqc::time_point time0;
  mpqc::time_point time1;
  time0 = mpqc::now(world_, accurate_time_);
  TArray result;

  BasisVector bs_array;
  std::shared_ptr<utility::TSPool<libint2::Engine>> engine_pool;
  std::shared_ptr<Screener> p_screener =
      std::make_shared<Screener>(Screener{});

  parse_two_body_three_center(formula, engine_pool, bs_array, p_screener);
  result = compute_integrals(this->world_, engine_pool, bs_array, p_screener);

  time1 = mpqc::now(world_, accurate_time_);
  time += mpqc::duration_in_s(time0, time1);

  utility::print_par(world_, "Computed Twobody Three Center Integral: ",
                     utility::to_string(formula.string()));
  double size = mpqc::detail::array_size(result);
  utility::print_par(world_, " Size: ", size, " GB");
  utility::print_par(world_, " Time: ", time, " s\n");

  return result;
}

template <typename Tile, typename Policy>
typename AOFactory<Tile, Policy>::TArray AOFactory<Tile, Policy>::compute4(
    const Formula& formula) {
  double time = 0.0;
  mpqc::time_point time0;
  mpqc::time_point time1;
  TArray result;

  if (formula.oper().has_option(Operator::Option::DensityFitting)) {
    // convert formula to df formula
    auto formula_strings = get_df_formula(formula);

    // compute integral
    auto left = compute(formula_strings[0]);
    auto center = compute(formula_strings[1]);
    auto right = compute(formula_strings[2]);

    time0 = mpqc::now(world_, accurate_time_);

    if (formula.notation() == Formula::Notation::Chemical) {
      result("i,j,k,l") = left("q,i,j") * center("q,p") * right("p,k,l");
    } else {
      result("i,j,k,l") = left("q,i,k") * center("q,p") * right("p,j,l");
    }

    time1 = mpqc::now(world_, accurate_time_);
    time += mpqc::duration_in_s(time0, time1);

    utility::print_par(
        world_, "Computed Twobody Four Center Density-Fitting Integral: ",
        utility::to_string(formula.string()));
    double size = mpqc::detail::array_size(result);
    utility::print_par(world_, " Size: ", size, " GB");
    utility::print_par(world_, " Time: ", time, " s\n");

  } else {
    time0 = mpqc::now(world_, accurate_time_);

    BasisVector bs_array;
    std::shared_ptr<Screener> p_screener =
        std::make_shared<Screener>(Screener{});

    std::shared_ptr<utility::TSPool<libint2::Engine>> engine_pool;
    parse_two_body_four_center(formula, engine_pool, bs_array, p_screener);
    result = compute_integrals(this->world_, engine_pool, bs_array, p_screener);

    // TODO handle permutation better
    if (formula.notation() == Formula::Notation::Physical) {
      result("i,j,k,l") = result("i,k,j,l");
    }

    time1 = mpqc::now(world_, accurate_time_);
    time += mpqc::duration_in_s(time0, time1);

    utility::print_par(world_, "Computed Twobody Four Center Integral: ",
                       utility::to_string(formula.string()));
    double size = mpqc::detail::array_size(result);
    utility::print_par(world_, " Size: ", size, " GB");
    utility::print_par(world_, " Time: ", time, " s\n");
  }
  return result;
}

extern template
class AOFactory<TA::TensorD, TA::SparsePolicy>;
//extern template
//class AOFactory<TA::TensorD, TA::DensePolicy>;

}  // namespace gaussian
}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_AO_FACTORY_H_
