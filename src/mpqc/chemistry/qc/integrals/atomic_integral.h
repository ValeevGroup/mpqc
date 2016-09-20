//
// Created by Chong Peng on 10/14/15.
//

#ifndef MPQC_ATOMIC_INTEGRAL_H
#define MPQC_ATOMIC_INTEGRAL_H

#include "../../../../../ta_routines/array_to_eigen.h"
#include "../../../../../ta_routines/sqrt_inv.h"
#include "../../../../../ta_routines/tile_convert.h"
#include "../../../../../utility/parallel_break_point.h"
#include "../../../../../utility/parallel_print.h"
#include "../../../../../utility/time.h"
#include "atomic_integral_base.h"
#include <madness/world/worldmem.h>
#include <mpqc/chemistry/qc/expression/permutation.h>
#include <mpqc/chemistry/qc/f12/f12_utility.h>
#include <mpqc/chemistry/qc/integrals/integrals.h>
#include <rapidjson/document.h>

namespace mpqc {
namespace integrals {

// TODO rename AtomicIntegral -> OperatorAOFactory
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

template <typename Tile, typename Policy>
class AtomicIntegral : public AtomicIntegralBase {
 public:
  using TArray = TA::DistArray<Tile, Policy>;

  /// Op is a function pointer that convert TA::Tensor to Tile
  //  using Op = Tile (*) (TA::TensorD&&);
  using Op = std::function<Tile(TA::TensorD&&)>;

  AtomicIntegral() = default;

  /**
   *  Constructor
   *  @param world madness::World object
   *  @param op op is a function that will take TA::TensorD as argument and
   *return Tile
   *  @param mol shared pointer to Molecule
   *  @param obs shared pointer to OrbitalBasisRegistry
   *  @param gtg_params  parameters used in computing f12 integrals
   *  @param in rapidjson Document object
   *
   *
   *
   *
   *  Options in Input
   *  @param AccurateTime, bool, control if use fence in timing, default false
   *  @param IterativeInvSqrt, bool, if use iterative approach to compute
   *inverse square root, defualt false
   */

  AtomicIntegral(madness::World& world, Op op,
                 const std::shared_ptr<molecule::Molecule>& mol,
                 const std::shared_ptr<basis::OrbitalBasisRegistry>& obs,
                 const std::vector<std::pair<double, double>>& gtg_params =
                     std::vector<std::pair<double, double>>(),
                 const rapidjson::Document& in = rapidjson::Document())
      : AtomicIntegralBase(world, mol, obs, gtg_params, in),
        ao_formula_registry_(),
        orbital_space_registry_(),
        op_(op) {
    if (in.IsObject()) {
      accurate_time_ =
          in.HasMember("AccurateTime") ? in["AccurateTime"].GetBool() : false;

      iterative_inv_sqrt_ = in.HasMember("IterativeInvSqrt")
                                ? in["IterativeInvSqrt"].GetBool()
                                : false;

    } else {
      accurate_time_ = false;
      iterative_inv_sqrt_ = false;
    }
    utility::print_par(world, "AccurateTime: ", accurate_time_, "\n");
    utility::print_par(world, "IterativeInvSqrt: ", iterative_inv_sqrt_, "\n");
  }

  /**
   * KeyVal Constructor
   *
   * It takes all the options from AtomicIntegralBase
   *
   * @param accurate_time, bool, control if use fence in timing, default false
   *
   * @return
   */

  AtomicIntegral(const KeyVal& kv)
      : AtomicIntegralBase(kv),
        ao_formula_registry_(),
        orbital_space_registry_() {
    accurate_time_ = kv.value("accurate_time", false);
    iterative_inv_sqrt_ = kv.value("iterative_inv_sqrt", false);

    /// Warning!!!!
    /// This is temporary workround
    /// For other Tile type, need a better way to set Op
    op_ = mpqc::ta_routines::TensorDPassThrough();
  }

  AtomicIntegral(AtomicIntegral&&) = default;
  AtomicIntegral& operator=(AtomicIntegral&&) = default;

  virtual ~AtomicIntegral() noexcept = default;

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
      Bvector const& bases,
      std::shared_ptr<Screener> p_screen =
          std::make_shared<integrals::Screener>(integrals::Screener{})) {
    auto result =
        mpqc::integrals::sparse_integrals(world, engine, bases, p_screen, op_);
    return result;
  }

  /// compute dense array
  template <typename U = Policy>
  TA::DistArray<Tile,
                typename std::enable_if<std::is_same<U, TA::DensePolicy>::value,
                                        TA::DensePolicy>::type>
  compute_integrals(
      madness::World& world, ShrPool<libint2::Engine>& engine,
      Bvector const& bases,
      std::shared_ptr<Screener> p_screen =
          std::make_shared<integrals::Screener>(integrals::Screener{})) {
    auto result =
        mpqc::integrals::dense_integrals(world, engine, bases, p_screen, op_);
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
typename AtomicIntegral<Tile, Policy>::TArray
AtomicIntegral<Tile, Policy>::compute(const std::wstring& formula_string) {
  auto formula = Formula(formula_string);
  return compute(formula);
}

template <typename Tile, typename Policy>
typename AtomicIntegral<Tile, Policy>::TArray
AtomicIntegral<Tile, Policy>::compute(const Formula& formula) {
  TArray result;

  // find if in registry
  auto iter = ao_formula_registry_.find(formula);

  if (iter != ao_formula_registry_.end()) {
    result = *(iter->second);
    utility::print_par(world_, "Retrieved AO Integral: ", utility::to_string(formula.string()));
    double size = utility::array_size(result);
    utility::print_par(world_, " Size: ", size, " GB\n");
    return result;
  } else {
    // find a permutation, currently, it won't store permutation in registry

    std::vector<Formula> permutes = permutations(formula);
    typename FormulaRegistry<TArray>::iterator find_permute;

    for (auto& permute : permutes) {
      find_permute = ao_formula_registry_.find(permute);
      if (find_permute != ao_formula_registry_.end()) {
        mpqc_time::t_point time0 = mpqc_time::now(world_, accurate_time_);

        // permute the array
        result(formula.to_ta_expression()) =
            (*(find_permute->second))(permute.to_ta_expression());

        mpqc_time::t_point time1 = mpqc_time::now(world_, accurate_time_);
        double time = mpqc_time::duration_in_s(time0, time1);

        utility::print_par(world_, "Permuted AO Integral: ",
                           utility::to_string(formula.string()), " From ",
                           utility::to_string(permute.string()));
        double size = utility::array_size(result);
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

    madness::print_meminfo(world_.rank(), "AOFactory: " + utility::to_string(formula.string()));

    return ao_formula_registry_.retrieve(formula);
  }
  // wait all process to obtain and insert result
  //        world_.gop.fence();
}

template <typename Tile, typename Policy>
typename AtomicIntegral<Tile, Policy>::TArray
AtomicIntegral<Tile, Policy>::compute2(const Formula& formula) {
  Bvector bs_array;
  double time = 0.0;
  mpqc_time::t_point time0;
  mpqc_time::t_point time1;
  TArray result;

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

      time0 = mpqc_time::now(world_, accurate_time_);

      result("i,j") = v("i,j") + t("i,j");

      time1 = mpqc_time::now(world_, accurate_time_);
      time += mpqc_time::duration_in_s(time0, time1);
    }
    // one body integral S, V, T...
    else {
      time0 = mpqc_time::now(world_, accurate_time_);

      std::shared_ptr<EnginePool<libint2::Engine>> engine_pool;
      parse_one_body(formula, engine_pool, bs_array);
      result = compute_integrals(this->world_, engine_pool, bs_array);

      time1 = mpqc_time::now(world_, accurate_time_);
      time += mpqc_time::duration_in_s(time0, time1);
    }
    utility::print_par(world_, "Computed One Body Integral: ",
                       utility::to_string(formula.string()));
    double size = utility::array_size(result);
    utility::print_par(world_, " Size: ", size, " GB");
    utility::print_par(world_, " Time: ", time, " s\n");
  }
  // use two body engine
  else if (formula.oper().is_twobody()) {
    time0 = mpqc_time::now(world_, accurate_time_);

    // compute inverse square root first in this case
    if (iterative_inv_sqrt_ &&
        formula.oper().has_option(Operator::Option::Inverse)) {
      auto inv_sqrt_formula = formula;
      inv_sqrt_formula.set_operator_option(
          {Operator::Option::InverseSquareRoot});

      result = this->compute(inv_sqrt_formula);

      time0 = mpqc_time::now(world_, accurate_time_);
      result("p,q") = result("p,r") * result("r,q");

      if (formula.oper().type() == Operator::Type::cGTG ||
          formula.oper().type() == Operator::Type::cGTGCoulomb) {
        result("i,j") = -result("i,j");
      }

    } else {
      // compute integral
      std::shared_ptr<EnginePool<libint2::Engine>> engine_pool;
      parse_two_body_two_center(formula, engine_pool, bs_array);
      result = compute_integrals(this->world_, engine_pool, bs_array);

      // inverse of integral
      if (formula.oper().has_option(Operator::Option::Inverse)) {
        if (formula.oper().type() == Operator::Type::cGTG ||
            formula.oper().type() == Operator::Type::cGTGCoulomb) {
          result("i,j") = -result("i,j");
        }

        //                utility::parallel_break_point(world_,0);
        //                std::cout << "Before Array To Eigen" << std::endl;
        //                std::cout << result << std::endl;
        MatrixD result_eig = array_ops::array_to_eigen(result);

        // compute cholesky decomposition
        auto llt_solver = Eig::LLT<MatrixD>(result_eig);

        // check success
        Eigen::ComputationInfo info = llt_solver.info();
        if (info == Eigen::ComputationInfo::Success) {
          MatrixD L = MatrixD(llt_solver.matrixL());
          MatrixD L_inv_eig = L.inverse();
          result_eig = L_inv_eig.transpose() * L_inv_eig;
        } else if (info == Eigen::ComputationInfo::NumericalIssue) {
          utility::print_par(
              world_,
              "!!!\nWarning!! NumericalIssue in Cholesky Decomposition\n!!!\n");
        } else if (info == Eigen::ComputationInfo::NoConvergence) {
          utility::print_par(
              world_,
              "!!!\nWarning!! NoConvergence in Cholesky Decomposition\n!!!\n");
        } else if (info == Eigen::ComputationInfo::InvalidInput) {
          utility::print_par(
              world_,
              "!!!\nWarning!! InvalidInput in Cholesky Decomposition\n!!!\n");
        }

        if (info != Eigen::ComputationInfo::Success) {
          utility::print_par(world_, "Using Eigen LU Decomposition Inverse!\n");

          Eigen::FullPivLU<MatrixD> lu(result_eig);

          TA_ASSERT(lu.isInvertible());

          result_eig = lu.inverse();
        }

        auto tr_result = result.trange().data()[0];
        result = array_ops::eigen_to_array<TA::TensorD>(
            result.get_world(), result_eig, tr_result, tr_result);

        if (formula.oper().type() == Operator::Type::cGTG ||
            formula.oper().type() == Operator::Type::cGTGCoulomb) {
          result("i,j") = -result("i,j");
        }
      }

      // inverse square root of integral
      if (formula.oper().has_option(Operator::Option::InverseSquareRoot)) {
        if (formula.oper().type() == Operator::Type::cGTG ||
            formula.oper().type() == Operator::Type::cGTGCoulomb) {
          result("i,j") = -result("i,j");
        }

        if (iterative_inv_sqrt_) {
          TArray tmp = array_ops::inverse_sqrt(result);
          result("i,j") = tmp("i,j");
        } else {
          auto result_eig = array_ops::array_to_eigen(result);
          MatrixD L_inv_eig =
              MatrixD(Eig::LLT<MatrixD>(result_eig).matrixL()).inverse();
          auto tr_result = result.trange().data()[0];
          result = array_ops::eigen_to_array<TA::TensorD>(
              result.get_world(), L_inv_eig, tr_result, tr_result);
        }

        if (formula.oper().type() == Operator::Type::cGTG ||
            formula.oper().type() == Operator::Type::cGTGCoulomb) {
          result("i,j") = -result("i,j");
        }
      }
    }

    time1 = mpqc_time::now(world_, accurate_time_);
    time += mpqc_time::duration_in_s(time0, time1);

    utility::print_par(world_, "Computed Twobody Two Center Integral: ",
                       utility::to_string(formula.string()));
    double size = utility::array_size(result);
    utility::print_par(world_, " Size: ", size, " GB");
    utility::print_par(world_, " Time: ", time, " s\n");
  }
  // compute JK, requires orbital space registry
  else if (formula.oper().is_jk()) {
    // density fitting case

    // find the density
    auto space_index = get_jk_orbital_space(formula.oper());
    auto& space = orbital_space_registry_->retrieve(space_index);

    auto obs = space.ao_key().name();
    if (formula.oper().has_option(Operator::Option::DensityFitting)) {
      auto three_center_formula = get_jk_df_formula(formula, obs);

      auto left = compute(three_center_formula[0]);
      auto center = compute(three_center_formula[1]);
      auto right = compute(three_center_formula[2]);

      time0 = mpqc_time::now(world_, accurate_time_);

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

      time1 = mpqc_time::now(world_, accurate_time_);
      time += mpqc_time::duration_in_s(time0, time1);
    }
    // four center case
    else {
      // find the density
      auto space_index = get_jk_orbital_space(formula.oper());
      auto& space = orbital_space_registry_->retrieve(space_index);
      auto obs = space.ao_key().name();
      // convert to ao formula
      auto four_center_formula = get_jk_formula(formula, obs);
      auto four_center = this->compute(four_center_formula);

      time0 = mpqc_time::now(world_, accurate_time_);

      if (formula.notation() == Formula::Notation::Chemical) {
        if (formula.oper().type() == Operator::Type::J) {
          result("rho,sigma") =
              four_center("rho,sigma,mu,nu") * (space("mu,i") * space("nu,i"));
        } else {
          result("rho,sigma") =
              four_center("rho,mu,sigma,nu") * (space("mu,i") * space("nu,i"));
        }
      } else {
        if (formula.oper().type() == Operator::Type::K) {
          result("rho,sigma") =
              four_center("rho,sigma,mu,nu") * (space("mu,i") * space("nu,i"));
        } else {
          result("rho,sigma") =
              four_center("rho,mu,sigma,nu") * (space("mu,i") * space("nu,i"));
        }
      }

      time1 = mpqc_time::now(world_, accurate_time_);
      time += mpqc_time::duration_in_s(time0, time1);
    }
    utility::print_par(world_, "Computed Coulumb/Exchange Integral: ",
                       utility::to_string(formula.string()));
    double size = utility::array_size(result);
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

    time0 = mpqc_time::now(world_, accurate_time_);

    result("i,j") = h("i,j") + 2 * j("i,j");

    time1 = mpqc_time::now(world_, accurate_time_);
    time += mpqc_time::duration_in_s(time0, time1);

    utility::print_par(world_, "Computed Coulumb/Exchange Integral: ",
                       utility::to_string(formula.string()));
    double size = utility::array_size(result);
    utility::print_par(world_, " Size: ", size, " GB");
    utility::print_par(world_, " Time: ", time, " s\n");
  }
  // compute Fock, requires orbital space registry
  else if (formula.oper().is_fock()) {
    auto formulas = get_fock_formula(formula);

    auto h = compute(formulas[0]);
    auto j = compute(formulas[1]);
    auto k = compute(formulas[2]);

    time0 = mpqc_time::now(world_, accurate_time_);
    // if closed shell
    if (formula.oper().type() == Operator::Type::Fock) {
      result("rho,sigma") =
          h("rho,sigma") + 2 * j("rho,sigma") - k("rho,sigma");
    }
    // else if spin orbital
    else {
      result("rho,sigma") = h("rho,sigma") + j("rho,sigma") - k("rho,sigma");
    }

    time1 = mpqc_time::now(world_, accurate_time_);

    time += mpqc_time::duration_in_s(time0, time1);

    utility::print_par(world_, "Computed Fock Integral: ",
                       utility::to_string(formula.string()));
    double size = utility::array_size(result);
    utility::print_par(world_, " Size: ", size, " GB");
    utility::print_par(world_, " Time: ", time, " s\n");
  }

  return result;
}

template <typename Tile, typename Policy>
typename AtomicIntegral<Tile, Policy>::TArray
AtomicIntegral<Tile, Policy>::compute3(const Formula& formula) {
  double time = 0.0;
  mpqc_time::t_point time0;
  mpqc_time::t_point time1;
  time0 = mpqc_time::now(world_, accurate_time_);
  TArray result;

  Bvector bs_array;
  std::shared_ptr<EnginePool<libint2::Engine>> engine_pool;
  std::shared_ptr<Screener> p_screener =
      std::make_shared<integrals::Screener>(integrals::Screener{});

  parse_two_body_three_center(formula, engine_pool, bs_array, p_screener);
  result = compute_integrals(this->world_, engine_pool, bs_array, p_screener);

  time1 = mpqc_time::now(world_, accurate_time_);
  time += mpqc_time::duration_in_s(time0, time1);

  utility::print_par(world_, "Computed Twobody Three Center Integral: ",
                     utility::to_string(formula.string()));
  double size = utility::array_size(result);
  utility::print_par(world_, " Size: ", size, " GB");
  utility::print_par(world_, " Time: ", time, " s\n");

  return result;
}

template <typename Tile, typename Policy>
typename AtomicIntegral<Tile, Policy>::TArray
AtomicIntegral<Tile, Policy>::compute4(const Formula& formula) {
  double time = 0.0;
  mpqc_time::t_point time0;
  mpqc_time::t_point time1;
  TArray result;

  if (formula.oper().has_option(Operator::Option::DensityFitting)) {
    // convert formula to df formula
    auto formula_strings = get_df_formula(formula);

    // compute integral
    auto left = compute(formula_strings[0]);
    auto center = compute(formula_strings[1]);
    auto right = compute(formula_strings[2]);

    time0 = mpqc_time::now(world_, accurate_time_);

    if (formula.notation() == Formula::Notation::Chemical) {
      result("i,j,k,l") = left("q,i,j") * center("q,p") * right("p,k,l");
    } else {
      result("i,j,k,l") = left("q,i,k") * center("q,p") * right("p,j,l");
    }

    time1 = mpqc_time::now(world_, accurate_time_);
    time += mpqc_time::duration_in_s(time0, time1);

    utility::print_par(
        world_, "Computed Twobody Four Center Density-Fitting Integral: ",
        utility::to_string(formula.string()));
    double size = utility::array_size(result);
    utility::print_par(world_, " Size: ", size, " GB");
    utility::print_par(world_, " Time: ", time, " s\n");

  } else {
    time0 = mpqc_time::now(world_, accurate_time_);

    Bvector bs_array;
    std::shared_ptr<Screener> p_screener =
        std::make_shared<integrals::Screener>(integrals::Screener{});

    std::shared_ptr<EnginePool<libint2::Engine>> engine_pool;
    parse_two_body_four_center(formula, engine_pool, bs_array, p_screener);
    result = compute_integrals(this->world_, engine_pool, bs_array, p_screener);

    // TODO handle permutation better
    if (formula.notation() == Formula::Notation::Physical) {
      result("i,j,k,l") = result("i,k,j,l");
    }

    time1 = mpqc_time::now(world_, accurate_time_);
    time += mpqc_time::duration_in_s(time0, time1);

    utility::print_par(world_, "Computed Twobody Four Center Integral: ",
                       utility::to_string(formula.string()));
    double size = utility::array_size(result);
    utility::print_par(world_, " Size: ", size, " GB");
    utility::print_par(world_, " Time: ", time, " s\n");
  }
  return result;
}
}
}

#endif  // MPQC_ATOMIC_INTEGRAL_H
