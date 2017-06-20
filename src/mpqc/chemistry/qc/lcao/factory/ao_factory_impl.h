//
// Created by Chong Peng on 2/2/17.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_AO_FACTORY_IMPL_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_AO_FACTORY_IMPL_H_

#include <regex>
#include <string>

#include "mpqc/math/groups/petite_list.h"
#include "mpqc/math/linalg/cholesky_inverse.h"

namespace mpqc {
namespace lcao {
namespace gaussian {

template <typename Tile, typename Policy>
AOFactory<Tile, Policy>::AOFactory(const KeyVal& kv)
    : Factory<TA::DistArray<Tile, Policy>, DirectArray<Tile, Policy>>(kv),
      gtg_params_() {
  ExEnv::out0() << "\nConstructing AOFactory: \n";
  std::string prefix = "";
  if (kv.exists("wfn_world") || kv.exists_class("wfn_world")) {
    prefix = "wfn_world:";
  }

  detail::set_oper<Tile>(op_);

  // if have auxilary basis
  if (kv.exists(prefix + "aux_basis")) {
    f12::GTGParams gtg_params;
    if (kv.exists(prefix + "f12_factor")) {
      auto f12_factor_str =
          kv.value<std::string>(prefix + "f12_factor", "stg-6g[1]");

      std::regex f12_factor_regex(
          "(stg|STG)\\-(\\d+)(g|G)\\[([0-9\\-eEdD\\.]+)\\]");
      std::smatch f12_factor_match;
      if (std::regex_search(f12_factor_str, f12_factor_match,
                            f12_factor_regex)) {
        auto n_gtg = std::stoi(f12_factor_match[2]);
        auto lengthscale = std::stod(f12_factor_match[4]);
        gtg_params = f12::GTGParams(lengthscale, n_gtg);
      } else
        throw InputError("invalid format for the f12_factor value", __FILE__,
                         __LINE__, "f12_factor");
    } else {
      if (kv.exists(prefix + "vir_basis")) {
        std::string basis_name =
            kv.value<std::string>(prefix + "vir_basis:name");
        gtg_params = f12::GTGParams(basis_name, 6);
      } else {
        std::string basis_name = kv.value<std::string>(prefix + "basis:name");
        gtg_params = f12::GTGParams(basis_name, 6);
      }
    }
    gtg_params_ = gtg_params.compute();
    ExEnv::out0() << indent
                  << "F12 Correlation Factor = " << gtg_params.exponent
                  << std::endl;
    ExEnv::out0() << indent << "NFunction = " << gtg_params.n_fit << std::endl;
    ExEnv::out0() << indent << "F12 Exponent Coefficient: \n";
    for (auto& pair : gtg_params_) {
      ExEnv::out0() << pair.first << " " << pair.second << std::endl;
    }
    ExEnv::out0() << std::endl;
  }
  // other initialization
  screen_ = kv.value<std::string>(prefix + "screen", "schwarz");
  screen_threshold_ = kv.value<double>(prefix + "screen_threshold", 1.0e-12);
  auto default_precision = std::numeric_limits<double>::epsilon();
  precision_ = kv.value<double>(prefix + "precision", default_precision);
  detail::integral_engine_precision = precision_;

  ExEnv::out0() << indent << "Screen = " << (screen_.empty() ? "none" : screen_)
                << "\n";
  if (!screen_.empty()) {
    ExEnv::out0() << indent << "ScreenThreshold = " << screen_threshold_
                  << "\n";
  }
  ExEnv::out0() << indent << "Precision = " << precision_ << "\n";
  iterative_inv_sqrt_ = kv.value<bool>(prefix + "iterative_inv_sqrt", false);
  ExEnv::out0() << indent << "Iterative inverse = "
                << (iterative_inv_sqrt_ ? "true" : "false") << "\n\n";
}

template <typename Tile, typename Policy>
typename AOFactory<Tile, Policy>::TArray AOFactory<Tile, Policy>::compute(
    const Formula& formula) {
  TA_USER_ASSERT(lcao::detail::if_all_ao(formula),
                 "AOFactory only accept AO index!\n");

  ExEnv::out0() << incindent;

  TArray result;

  auto& world = this->world();

  // find if in registry
  auto iter = this->registry_.find(formula);

  if (iter != this->registry_.end()) {
    result = iter->second;
    if(this->verbose_){
      ExEnv::out0() << indent;
      ExEnv::out0() << "Retrieved AO Integral: "
                    << utility::to_string(formula.string());
      double size = mpqc::detail::array_size(result);
      ExEnv::out0() << " Size: " << size << " GB\n";
    }
  } else {
    // find a permutation, currently, it won't store permutation in registry

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

        if(this->verbose_){
          ExEnv::out0() << indent;
          ExEnv::out0() << "Permuted AO Integral: "
                        << utility::to_string(formula.string()) << " From "
                        << utility::to_string(permute.string());
          double size = mpqc::detail::array_size(result);
          ExEnv::out0() << " Size: " << size << " GB"
                        << " Time: " << time << " s\n";
        }

        // store current array and delete old one
        this->registry_.insert(formula, result);
        this->registry_.purge_formula(permute);
      }
    }

    // if formula not found
    if (!result.is_initialized()) {
      // compute formula
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
        world.rank(), "AOFactory: " + utility::to_string(formula.string()));
  }
  ExEnv::out0() << decindent;
  return result;
}

template <typename Tile, typename Policy>
typename AOFactory<Tile, Policy>::TArray AOFactory<Tile, Policy>::compute2(
    const Formula& formula) {
  BasisVector bs_array;
  double time = 0.0;
  mpqc::time_point time0;
  mpqc::time_point time1;
  TArray result;
  auto& world = this->world();

  // get the inverse square root instead
  if (iterative_inv_sqrt_ && formula.has_option(Formula::Option::Inverse)) {
    auto inv_sqrt_formula = formula;
    inv_sqrt_formula.clear_option();
    inv_sqrt_formula.add_option(Formula::Option::InverseSquareRoot);

    result = this->compute(inv_sqrt_formula);

    time0 = mpqc::now(world, this->accurate_time_);
    result("p,q") = result("p,r") * result("r,q");

    if (formula.oper().type() == Operator::Type::cGTG ||
        formula.oper().type() == Operator::Type::cGTGCoulomb) {
      result("i,j") = -result("i,j");
    }

    if(this->verbose_){
      ExEnv::out0() << indent;
      ExEnv::out0() << "Computed Inverse of Integral: "
                    << utility::to_string(formula.string());
      double size = mpqc::detail::array_size(result);

      time1 = mpqc::now(world, this->accurate_time_);
      time += mpqc::duration_in_s(time0, time1);
      ExEnv::out0() << " Size: " << size << " GB"
                    << " Time: " << time << " s\n";
    }
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

        time0 = mpqc::now(world, this->accurate_time_);

        result("i,j") = v("i,j") + t("i,j");

        time1 = mpqc::now(world, this->accurate_time_);
        time += mpqc::duration_in_s(time0, time1);
      }
      // one body integral S, V, T...
      else {
        time0 = mpqc::now(world, this->accurate_time_);

        std::shared_ptr<utility::TSPool<libint2::Engine>> engine_pool;
        parse_one_body(formula, engine_pool, bs_array);
        result = compute_integrals(world, engine_pool, bs_array);

        time1 = mpqc::now(world, this->accurate_time_);
        time += mpqc::duration_in_s(time0, time1);
      }

      if(this->verbose_){
        ExEnv::out0() << indent;
        ExEnv::out0() << "Computed One Body Integral: "
                      << utility::to_string(formula.string());
        double size = mpqc::detail::array_size(result);
        ExEnv::out0() << " Size: " << size << " GB"
                      << " Time: " << time << " s\n";
      }
    }
    // use two body engine
    else if (formula.oper().is_twobody()) {
      time0 = mpqc::now(world, this->accurate_time_);

      // compute integral
      std::shared_ptr<utility::TSPool<libint2::Engine>> engine_pool;
      parse_two_body_two_center(formula, engine_pool, bs_array);
      result = compute_integrals(world, engine_pool, bs_array);

      time1 = mpqc::now(world, this->accurate_time_);
      time += mpqc::duration_in_s(time0, time1);

      if(this->verbose_){
        ExEnv::out0() << indent;
        ExEnv::out0() << "Computed Twobody Two Center Integral: "
                      << utility::to_string(formula.string());
        double size = mpqc::detail::array_size(result);
        ExEnv::out0() << " Size: " << size << " GB"
                      << " Time: " << time << " s\n";
      }
    }
    // compute JK, requires orbital space registry
    else if (formula.oper().is_jk()) {
      // density fitting case

      // find the density
      auto space_index = detail::get_jk_orbital_space(formula.oper());
      auto& space = this->orbital_registry().retrieve(space_index);

      auto obs = space.ao_index().name();
      if (formula.has_option(Formula::Option::DensityFitting)) {
        auto three_center_formula = detail::get_jk_df_formula(formula, obs);

        auto left = compute_direct(three_center_formula[0]);
        auto center = compute(three_center_formula[1]);
        auto right = compute_direct(three_center_formula[2]);

        time0 = mpqc::now(world, this->accurate_time_);

        // J case
        if (formula.oper().type() == Operator::Type::J) {
          result("i,j") = space("k,a") * space("l,a") * right("Q,k,l") *
                          center("K,Q") * left("K,i,j");
        }
        // K case
        else {
          result("i,j") = (left("K,i,k") * space("k,a")) * center("K,Q") *
                          (right("Q,j,l") * space("l,a"));
        }

        time1 = mpqc::now(world, this->accurate_time_);
        time += mpqc::duration_in_s(time0, time1);
      }
      // four center case
      else {
        // find the density
        auto space_index = detail::get_jk_orbital_space(formula.oper());
        auto& space = this->orbital_registry().retrieve(space_index);
        auto obs = space.ao_index().name();
        // convert to ao formula
        auto four_center_formula = detail::get_jk_formula(formula, obs);
        auto four_center = this->compute(four_center_formula);

        time0 = mpqc::now(world, this->accurate_time_);

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

        time1 = mpqc::now(world, this->accurate_time_);
        time += mpqc::duration_in_s(time0, time1);
      }

      if(this->verbose_){
        ExEnv::out0() << indent;
        ExEnv::out0() << "Computed Coulumb/Exchange Integral: "
                      << utility::to_string(formula.string());
        double size = mpqc::detail::array_size(result);
        ExEnv::out0() << " Size: " << size << " GB"
                      << " Time: " << time << " s\n";
      }
    }
    // hJ = H + J
    else if (formula.oper().type() == Operator::Type::hJ) {
      auto h_formula = formula;
      h_formula.set_operator(Operator(L"H"));

      auto j_formula = formula;
      j_formula.set_operator_type(Operator::Type::J);

      auto h = this->compute(h_formula);
      auto j = this->compute(j_formula);

      time0 = mpqc::now(world, this->accurate_time_);

      result("i,j") = h("i,j") + 2 * j("i,j");

      time1 = mpqc::now(world, this->accurate_time_);
      time += mpqc::duration_in_s(time0, time1);

      if(this->verbose_){
        ExEnv::out0() << indent;
        ExEnv::out0() << "Computed Coulumb/Exchange Integral: "
                      << utility::to_string(formula.string());
        double size = mpqc::detail::array_size(result);
        ExEnv::out0() << " Size: " << size << " GB"
                      << " Time: " << time << " s\n";
      }
    }
    // compute Fock, requires orbital space registry
    else if (formula.oper().is_fock()) {
      auto formulas = detail::get_fock_formula(formula);

      auto h = compute(formulas[0]);
      auto j = compute(formulas[1]);
      auto k = compute(formulas[2]);

      time0 = mpqc::now(world, this->accurate_time_);
      // if closed shell
      if (formula.oper().type() == Operator::Type::Fock) {
        result("rho,sigma") =
            h("rho,sigma") + 2 * j("rho,sigma") - k("rho,sigma");
      }
      // else if spin orbital
      else {
        result("rho,sigma") = h("rho,sigma") + j("rho,sigma") - k("rho,sigma");
      }

      time1 = mpqc::now(world, this->accurate_time_);

      time += mpqc::duration_in_s(time0, time1);

      if(this->verbose_){

        ExEnv::out0() << indent;
        ExEnv::out0() << "Computed Fock Integral: "
                      << utility::to_string(formula.string());
        double size = mpqc::detail::array_size(result);
        ExEnv::out0() << " Size: " << size << " GB"
                      << " Time: " << time << " s\n";
      }
    }
  }

  // check for other options

  // compute inverse
  if (!iterative_inv_sqrt_ && formula.has_option(Formula::Option::Inverse)) {
    time0 = mpqc::now(world, this->accurate_time_);

    if (formula.oper().type() == Operator::Type::cGTG ||
        formula.oper().type() == Operator::Type::cGTGCoulomb) {
      result("i,j") = -result("i,j");
    }

    auto tmp = array_ops::eigen_inverse(result);
    result("i,j") = tmp("i,j");

    if (formula.oper().type() == Operator::Type::cGTG ||
        formula.oper().type() == Operator::Type::cGTGCoulomb) {
      result("i,j") = -result("i,j");
    }
    time1 = mpqc::now(world, this->accurate_time_);
    auto inv_time = mpqc::duration_in_s(time0, time1);
    if(this->verbose_){
      ExEnv::out0() << indent;
      ExEnv::out0() << "Inverse Time: " << inv_time << " s\n";
    }
  }

  // inverse square root of integral
  if (formula.has_option(Formula::Option::InverseSquareRoot)) {
    time0 = mpqc::now(world, this->accurate_time_);
    if (formula.oper().type() == Operator::Type::cGTG ||
        formula.oper().type() == Operator::Type::cGTGCoulomb) {
      result("i,j") = -result("i,j");
    }

    if (iterative_inv_sqrt_) {
      TArray tmp = array_ops::inverse_sqrt(result);
      result("i,j") = tmp("i,j");
    } else {
      auto result_eig = array_ops::array_to_eigen(result);

      Eigen::SelfAdjointEigenSolver<decltype(result_eig)> es(result_eig);
      decltype(result_eig) inv_eig = es.operatorInverseSqrt();

      auto tr_result = result.trange().data()[0];
      result = array_ops::eigen_to_array<Tile, Policy>(result.world(), inv_eig,
                                                       tr_result, tr_result);
    }

    if (formula.oper().type() == Operator::Type::cGTG ||
        formula.oper().type() == Operator::Type::cGTGCoulomb) {
      result("i,j") = -result("i,j");
    }
    time1 = mpqc::now(world, this->accurate_time_);
    auto inv_sqrt_time = mpqc::duration_in_s(time0, time1);
    if(this->verbose_){
      ExEnv::out0() << indent;
      ExEnv::out0() << "Inverse Square Root Time: " << inv_sqrt_time << " s\n";
    }
  }

  return result;
}
template <typename Tile, typename Policy>
typename AOFactory<Tile, Policy>::TArray
AOFactory<Tile, Policy>::compute_cadf_coeffs(const Formula& formula) {
  double time = 0.0;
  mpqc::time_point time0;
  mpqc::time_point time1;
  auto& world = this->world();
  time0 = mpqc::now(world, this->accurate_time_);

  auto Xindex = formula.bra_indices();
  auto mu_nu_index = formula.ket_indices();

  const auto& basis_registry = *this->basis_registry();

  auto obs = detail::index_to_basis(basis_registry, mu_nu_index[0]);
  auto dfbs = detail::index_to_basis(basis_registry, Xindex[0]);

  auto C = cadf_fitting_coefficients<Tile, Policy>(world, *obs, *dfbs);

  time1 = mpqc::now(world, this->accurate_time_);
  time += mpqc::duration_in_s(time0, time1);

  if(this->verbose_){
    ExEnv::out0() << indent;
    ExEnv::out0() << "Computed CADF fitting Coefficients: "
                  << utility::to_string(formula.string());
    double size = mpqc::detail::array_size(C);
    ExEnv::out0() << " Size: " << size << " GB"
                  << " Time: " << time << " s\n";
  }
  return C;
}

template <typename Tile, typename Policy>
typename AOFactory<Tile, Policy>::TArray AOFactory<Tile, Policy>::compute3(
    const Formula& formula) {


  // no option supported now
  TA_ASSERT(formula.options().empty());

  double time = 0.0;
  mpqc::time_point time0;
  mpqc::time_point time1;
  auto& world = this->world();
  time0 = mpqc::now(world, this->accurate_time_);
  TArray result;

  BasisVector bs_array;
  std::shared_ptr<utility::TSPool<libint2::Engine>> engine_pool;
  std::shared_ptr<Screener> p_screener = std::make_shared<Screener>(Screener{});

  if (formula.oper().type() == Operator::Type::Cadf) {
    return compute_cadf_coeffs(formula);
  }

  parse_two_body_three_center(formula, engine_pool, bs_array, p_screener);
  result = compute_integrals(this->world(), engine_pool, bs_array, p_screener);

  time1 = mpqc::now(world, this->accurate_time_);
  time += mpqc::duration_in_s(time0, time1);

  if(this->verbose_){
    ExEnv::out0() << indent;
    ExEnv::out0() << "Computed Twobody Three Center Integral: "
                  << utility::to_string(formula.string());
    double size = mpqc::detail::array_size(result);
    ExEnv::out0() << " Size: " << size << " GB"
                  << " Time: " << time << " s\n";
  }

  return result;
}

template <typename Tile, typename Policy>
typename AOFactory<Tile, Policy>::TArray AOFactory<Tile, Policy>::compute4(
    const Formula& formula) {
  double time = 0.0;
  mpqc::time_point time0;
  mpqc::time_point time1;
  auto& world = this->world();
  TArray result;

  if (formula.has_option(Formula::Option::DensityFitting)) {
    // convert formula to df formula
    auto formula_strings = detail::get_df_formula(formula);

    // compute integral
    auto left = compute(formula_strings[0]);
    auto right = compute(formula_strings[1]);

    time0 = mpqc::now(world, this->accurate_time_);

    std::string result_str = formula.notation() == Formula::Notation::Chemical
                                 ? "i,j,k,l"
                                 : "i,k,j,l";

    if (formula.oper().type() == Operator::Type::cGTG ||
        formula.oper().type() == Operator::Type::cGTGCoulomb)
    {
      result(result_str) = left("q,i,j") * (-right("q,k,l"));
    } else {
      result(result_str) = left("q,i,j") * right("q,k,l");
    }

    time1 = mpqc::now(world, this->accurate_time_);
    time += mpqc::duration_in_s(time0, time1);

    if(this->verbose_){
      ExEnv::out0() << indent;
      ExEnv::out0() << "Computed Twobody Four Center Density-Fitting Integral: "
                    << utility::to_string(formula.string());
      double size = mpqc::detail::array_size(result);
      ExEnv::out0() << " Size: " << size << " GB"
                    << " Time: " << time << " s\n";
    }

  } else {
    time0 = mpqc::now(world, this->accurate_time_);

    BasisVector bs_array;
    std::shared_ptr<Screener> p_screener =
        std::make_shared<Screener>(Screener{});

    std::shared_ptr<utility::TSPool<libint2::Engine>> engine_pool;
    parse_two_body_four_center(formula, engine_pool, bs_array, p_screener);
    result = compute_integrals(world, engine_pool, bs_array, p_screener);

    // TODO handle permutation better
    if (formula.notation() == Formula::Notation::Physical) {
      result("i,j,k,l") = result("i,k,j,l");
    }

    time1 = mpqc::now(world, this->accurate_time_);
    time += mpqc::duration_in_s(time0, time1);

    if(this->verbose_){
      ExEnv::out0() << indent;
      ExEnv::out0() << "Computed Twobody Four Center Integral: "
                    << utility::to_string(formula.string());
      double size = mpqc::detail::array_size(result);
      ExEnv::out0() << " Size: " << size << " GB"
                    << " Time: " << time << " s\n";
    }
  }
  return result;
}

///
///
///  Functions to compute direct integral
///
///

template <typename Tile, typename Policy>
typename AOFactory<Tile, Policy>::DirectTArray
AOFactory<Tile, Policy>::compute_direct(const Formula& formula) {
  TA_USER_ASSERT(lcao::detail::if_all_ao(formula),
                 "AOFactory only accept AO index!\n");

  ExEnv::out0() << incindent;

  DirectTArray result;

  // find if in registry
  auto iter = this->direct_registry_.find(formula);

  if (iter != this->direct_registry_.end()) {
    result = iter->second;

    if(this->verbose_){
      ExEnv::out0() << indent << "Retrieved Direct AO Integral: "
                    << utility::to_string(formula.string());
      double size = mpqc::detail::array_size(result.array());
      ExEnv::out0() << " Size: " << size << " GB\n";
    }

  } else {
    // compute formula
    if (formula.rank() == 2) {
      result = compute_direct2(formula);
      this->direct_registry_.insert(formula, result);
    } else if (formula.rank() == 3) {
      result = compute_direct3(formula);
      this->direct_registry_.insert(formula, result);
    } else if (formula.rank() == 4) {
      result = compute_direct4(formula);
      this->direct_registry_.insert(formula, result);
    }
  }

  ExEnv::out0() << decindent;
  return result;
}

template <typename Tile, typename Policy>
typename AOFactory<Tile, Policy>::DirectTArray
AOFactory<Tile, Policy>::compute_direct2(const Formula& formula) {
  BasisVector bs_array;
  double time = 0.0;
  std::shared_ptr<utility::TSPool<libint2::Engine>> engine_pool;
  mpqc::time_point time0;
  mpqc::time_point time1;
  auto& world = this->world();

  DirectTArray result;

  if (formula.oper().is_onebody()) {
    time0 = mpqc::now(world, this->accurate_time_);

    parse_one_body(formula, engine_pool, bs_array);
    result = compute_direct_integrals(world, engine_pool, bs_array);

    time1 = mpqc::now(world, this->accurate_time_);
    time += mpqc::duration_in_s(time0, time1);

    if(this->verbose_){
      ExEnv::out0() << indent << "Computed Direct One Body Integral: "
                    << utility::to_string(formula.string());
      double size = mpqc::detail::array_size(result.array());
      ExEnv::out0() << " Size: " << size << " GB"
                    << " Time: " << time << " s\n";
    }

  } else if (formula.oper().is_twobody()) {
    time0 = mpqc::now(world, this->accurate_time_);

    parse_two_body_two_center(formula, engine_pool, bs_array);
    result = compute_direct_integrals(world, engine_pool, bs_array);

    time1 = mpqc::now(world, this->accurate_time_);
    time += mpqc::duration_in_s(time0, time1);

    if(this->verbose_){
      ExEnv::out0() << indent << "Computed Direct Twobody Two Center Integral: "
                    << utility::to_string(formula.string());
      double size = mpqc::detail::array_size(result.array());
      ExEnv::out0() << " Size: " << size << " GB Time: " << time << " s\n";
    }
  } else {
    throw ProgrammingError("Unsupported Operator in DirectAOFactory!!\n",
                           __FILE__, __LINE__);
  }

  madness::print_meminfo(
      world.rank(), utility::wconcat("DirectAOFactory:", formula.string()));
  return result;
}

template <typename Tile, typename Policy>
typename AOFactory<Tile, Policy>::DirectTArray
AOFactory<Tile, Policy>::compute_direct3(const Formula& formula) {
  double time = 0.0;
  mpqc::time_point time0;
  mpqc::time_point time1;
  auto& world = this->world();
  time0 = mpqc::now(world, this->accurate_time_);
  DirectTArray result;

  BasisVector bs_array;
  std::shared_ptr<utility::TSPool<libint2::Engine>> engine_pool;
  std::shared_ptr<Screener> p_screener = std::make_shared<Screener>(Screener{});

  parse_two_body_three_center(formula, engine_pool, bs_array, p_screener);
  result = compute_direct_integrals(world, engine_pool, bs_array, p_screener);

  time1 = mpqc::now(world, this->accurate_time_);
  time += mpqc::duration_in_s(time0, time1);

  if(this->verbose_){
    ExEnv::out0() << indent << "Computed Direct Twobody Three Center Integral: "
                  << utility::to_string(formula.string());
    double size = mpqc::detail::array_size(result.array());
    ExEnv::out0() << " Size: " << size << " GB Time: " << time << " s\n";
  }
  madness::print_meminfo(
      world.rank(), utility::wconcat("DirectAOFactory:", formula.string()));

  return result;
}

template <typename Tile, typename Policy>
typename AOFactory<Tile, Policy>::DirectTArray
AOFactory<Tile, Policy>::compute_direct4(const Formula& formula) {
  if (formula.notation() != Formula::Notation::Chemical) {
    throw ProgrammingError(
        "Direct AO Integral Only Support Chemical Notation! \n", __FILE__,
        __LINE__);
  }

  double time = 0.0;
  mpqc::time_point time0;
  mpqc::time_point time1;
  auto& world = this->world();
  DirectTArray result;

  time0 = mpqc::now(world, this->accurate_time_);

  BasisVector bs_array;
  std::shared_ptr<Screener> p_screener = std::make_shared<Screener>(Screener{});
  std::shared_ptr<utility::TSPool<libint2::Engine>> engine_pool;

  parse_two_body_four_center(formula, engine_pool, bs_array, p_screener);
  auto plist = math::PetiteList::make(formula.symmetry());

  result =
      compute_direct_integrals(world, engine_pool, bs_array, p_screener, plist);

  time1 = mpqc::now(world, this->accurate_time_);
  time += mpqc::duration_in_s(time0, time1);

  if(this->verbose_){
    ExEnv::out0() << indent << "Computed Direct Twobody Four Center Integral: "
                  << utility::to_string(formula.string());
    double size = mpqc::detail::array_size(result.array());
    ExEnv::out0() << " Size: " << size << " GB Time: " << time << " s\n";
  }
  madness::print_meminfo(world.rank(),
                         utility::wconcat("AOFactory:", formula.string()));
  return result;
}

template <typename Tile, typename Policy>
void AOFactory<Tile, Policy>::parse_one_body(
    const Formula& formula,
    std::shared_ptr<utility::TSPool<libint2::Engine>>& engine_pool,
    BasisVector& bases) {
  auto bra_indices = formula.bra_indices();
  auto ket_indices = formula.ket_indices();

  TA_ASSERT(bra_indices.size() == 1);
  TA_ASSERT(ket_indices.size() == 1);

  auto bra_index = bra_indices[0];
  auto ket_index = ket_indices[0];

  TA_ASSERT(bra_index.is_ao());
  TA_ASSERT(ket_index.is_ao());

  const auto& basis_registry = *this->basis_registry();

  auto bra_basis = detail::index_to_basis(basis_registry, bra_index);
  auto ket_basis = detail::index_to_basis(basis_registry, ket_index);

  TA_ASSERT(bra_basis != nullptr);
  TA_ASSERT(ket_basis != nullptr);

  bases = BasisVector{{*bra_basis, *ket_basis}};

  auto oper_type = formula.oper().type();
  engine_pool = make_engine_pool(
      detail::to_libint2_operator(oper_type),
      utility::make_array_of_refs(*bra_basis, *ket_basis), libint2::BraKet::x_x,
      detail::to_libint2_operator_params(oper_type, *this->atoms(),
                                         this->gtg_params_));
}

template <typename Tile, typename Policy>
void AOFactory<Tile, Policy>::parse_two_body_two_center(
    const Formula& formula,
    std::shared_ptr<utility::TSPool<libint2::Engine>>& engine_pool,
    BasisVector& bases) {
  TA_USER_ASSERT(formula.notation() == Formula::Notation::Chemical,
                 "Two Body Two Center Integral Must Use Chemical Notation");

  auto bra_indexs = formula.bra_indices();
  auto ket_indexs = formula.ket_indices();

  TA_ASSERT(bra_indexs.size() == 1);
  TA_ASSERT(ket_indexs.size() == 1);

  auto bra_index0 = bra_indexs[0];
  auto ket_index0 = ket_indexs[0];

  TA_ASSERT(ket_index0.is_ao());
  TA_ASSERT(bra_index0.is_ao());

  const auto& basis_registry = *this->basis_registry();

  auto bra_basis0 = detail::index_to_basis(basis_registry, bra_index0);
  auto ket_basis0 = detail::index_to_basis(basis_registry, ket_index0);

  TA_ASSERT(bra_basis0 != nullptr);
  TA_ASSERT(ket_basis0 != nullptr);

  bases = BasisVector{{*bra_basis0, *ket_basis0}};

  auto oper_type = formula.oper().type();
  engine_pool =
      make_engine_pool(detail::to_libint2_operator(oper_type),
                       utility::make_array_of_refs(*bra_basis0, *ket_basis0),
                       libint2::BraKet::xs_xs,
                       detail::to_libint2_operator_params(
                           oper_type, *this->atoms(), this->gtg_params_));
}

template <typename Tile, typename Policy>
void AOFactory<Tile, Policy>::parse_two_body_three_center(
    const Formula& formula,
    std::shared_ptr<utility::TSPool<libint2::Engine>>& engine_pool,
    BasisVector& bases, std::shared_ptr<Screener>& p_screener) {
  TA_USER_ASSERT(formula.notation() == Formula::Notation::Chemical,
                 "Three Center Integral Must Use Chemical Notation");

  auto bra_indexs = formula.bra_indices();
  auto ket_indexs = formula.ket_indices();

  TA_ASSERT(bra_indexs.size() == 1);
  TA_ASSERT(ket_indexs.size() == 2);

  TA_ASSERT(bra_indexs[0].is_ao());
  TA_ASSERT(ket_indexs[0].is_ao());
  TA_ASSERT(ket_indexs[1].is_ao());

  const auto& basis_registry = *this->basis_registry();

  auto bra_basis0 = detail::index_to_basis(basis_registry, bra_indexs[0]);
  auto ket_basis0 = detail::index_to_basis(basis_registry, ket_indexs[0]);
  auto ket_basis1 = detail::index_to_basis(basis_registry, ket_indexs[1]);

  TA_ASSERT(bra_basis0 != nullptr);
  TA_ASSERT(ket_basis0 != nullptr);
  TA_ASSERT(ket_basis1 != nullptr);

  bases = BasisVector{{*bra_basis0, *ket_basis0, *ket_basis1}};

  auto oper_type = formula.oper().type();
  engine_pool = make_engine_pool(
      detail::to_libint2_operator(oper_type),
      utility::make_array_of_refs(*bra_basis0, *ket_basis0, *ket_basis1),
      libint2::BraKet::xs_xx,
      detail::to_libint2_operator_params(oper_type, *this->atoms(),
                                         this->gtg_params_));

  if (!screen_.empty()) {
    /// make another engine to screener!!!

    auto screen_engine_pool = make_engine_pool(
        detail::to_libint2_operator(oper_type),
        utility::make_array_of_refs(*bra_basis0, *ket_basis0, *ket_basis1),
        libint2::BraKet::xx_xx,
        detail::to_libint2_operator_params(oper_type, *this->atoms(),
                                           this->gtg_params_));

    p_screener = detail::make_screener(this->world(), screen_engine_pool, bases,
                                       this->screen_, this->screen_threshold_);
  }
}

template <typename Tile, typename Policy>
void AOFactory<Tile, Policy>::parse_two_body_four_center(
    const Formula& formula,
    std::shared_ptr<utility::TSPool<libint2::Engine>>& engine_pool,
    BasisVector& bases, std::shared_ptr<Screener>& p_screener) {
  auto bra_indexs = formula.bra_indices();
  auto ket_indexs = formula.ket_indices();

  TA_ASSERT(bra_indexs.size() == 2);
  TA_ASSERT(ket_indexs.size() == 2);

  TA_ASSERT(bra_indexs[0].is_ao());
  TA_ASSERT(ket_indexs[0].is_ao());
  TA_ASSERT(bra_indexs[1].is_ao());
  TA_ASSERT(ket_indexs[1].is_ao());

  const auto& basis_registry = *this->basis_registry();
  auto bra_basis0 = detail::index_to_basis(basis_registry, bra_indexs[0]);
  auto ket_basis0 = detail::index_to_basis(basis_registry, ket_indexs[0]);
  auto bra_basis1 = detail::index_to_basis(basis_registry, bra_indexs[1]);
  auto ket_basis1 = detail::index_to_basis(basis_registry, ket_indexs[1]);

  TA_ASSERT(bra_basis0 != nullptr);
  TA_ASSERT(ket_basis0 != nullptr);
  TA_ASSERT(bra_basis1 != nullptr);
  TA_ASSERT(ket_basis1 != nullptr);

  if (formula.notation() == Formula::Notation::Chemical) {
    bases = {{*bra_basis0, *bra_basis1, *ket_basis0, *ket_basis1}};
  } else {
    bases = {{*bra_basis0, *ket_basis0, *bra_basis1, *ket_basis1}};
  }

  auto oper_type = formula.oper().type();
  engine_pool = make_engine_pool(
      detail::to_libint2_operator(oper_type),
      utility::make_array_of_refs(bases[0], bases[1], bases[2], bases[3]),
      libint2::BraKet::xx_xx,
      detail::to_libint2_operator_params(oper_type, *this->atoms(),
                                         this->gtg_params_));

  p_screener = detail::make_screener(this->world(), engine_pool, bases,
                                     this->screen_, this->screen_threshold_);
}

}  // namespace gaussian

}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_AO_FACTORY_IMPL_H_
