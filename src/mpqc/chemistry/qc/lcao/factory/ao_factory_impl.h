//
// Created by Chong Peng on 2/2/17.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_AO_FACTORY_IMPL_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_AO_FACTORY_IMPL_H_

namespace mpqc {
namespace lcao {
namespace gaussian {

template <typename Tile, typename Policy>
AOFactory<Tile, Policy>::AOFactory(const KeyVal& kv)
    : world_(*kv.value<madness::World*>("$:world")),
      orbital_basis_registry_(),
      molecule_(),
      gtg_params_(),
      ao_formula_registry_(),
      orbital_space_registry_() {
  std::string prefix = "";
  if (kv.exists("wfn_world") || kv.exists_class("wfn_world")) {
    prefix = "wfn_world:";
  }

  accurate_time_ = kv.value<bool>(prefix + "accurate_time", false);
  iterative_inv_sqrt_ = kv.value<bool>(prefix + "iterative_inv_sqrt", false);

  set_oper(Tile());

  /// Basis will come from wfn_world
  //  orbital_basis_registry_ =
  //  std::make_shared<basis::OrbitalBasisRegistry>(basis::OrbitalBasisRegistry(kv));
  molecule_ = kv.class_ptr<Molecule>(prefix + "molecule");

  // if have auxilary basis
  if (kv.exists(prefix + "aux_basis")) {
    int n_function = kv.value<int>(prefix + "corr_functions", 6);
    double corr_param = kv.value<double>(prefix + "corr_param", 0);
    f12::GTGParams gtg_params;
    if (corr_param != 0) {
      gtg_params = f12::GTGParams(corr_param, n_function);
    } else {
      if (kv.exists(prefix + "vir_basis")) {
        std::string basis_name =
            kv.value<std::string>(prefix + "vir_basis:name");
        gtg_params = f12::GTGParams(basis_name, n_function);
      } else {
        std::string basis_name = kv.value<std::string>(prefix + "basis:name");
        gtg_params = f12::GTGParams(basis_name, n_function);
      }
    }
    gtg_params_ = gtg_params.compute();
    if (world().rank() == 0) {
      std::cout << "F12 Correlation Factor: " << gtg_params.exponent
                << std::endl;
      std::cout << "NFunction: " << gtg_params.n_fit << std::endl;
      std::cout << "F12 Exponent Coefficient" << std::endl;
      for (auto& pair : gtg_params_) {
        std::cout << pair.first << " " << pair.second << std::endl;
      }
      std::cout << std::endl;
    }
  }
  // other initialization
  screen_ = kv.value<std::string>(prefix + "screen", "");
  screen_threshold_ = kv.value<double>(prefix + "threshold", 1.0e-10);
  auto default_precision = std::numeric_limits<double>::epsilon();
  precision_ = kv.value<double>(prefix + "precision", default_precision);
  detail::integral_engine_precision = precision_;

  utility::print_par(world_, "Screen: ", screen_, "\n");
  if (!screen_.empty()) {
    utility::print_par(world_, "Threshold: ", screen_threshold_, "\n");
  }
  utility::print_par(world_, "Precision: ", precision_, "\n");
  utility::print_par(world_, "\n");
}

template <typename Tile, typename Policy>
typename AOFactory<Tile, Policy>::TArray AOFactory<Tile, Policy>::compute(
    const std::wstring& formula_string) {
  auto formula = Formula(formula_string);
  return compute(formula);
}

template <typename Tile, typename Policy>
typename AOFactory<Tile, Policy>::TArray AOFactory<Tile, Policy>::compute(
    const Formula& formula) {
  ExEnv::out0() << incindent;

  TArray result;

  // find if in registry
  auto iter = ao_formula_registry_.find(formula);

  if (iter != ao_formula_registry_.end()) {
    result = iter->second;
    ExEnv::out0() << indent;
    ExEnv::out0() << "Retrieved AO Integral: "
                  << utility::to_string(formula.string());
    double size = mpqc::detail::array_size(result);
    ExEnv::out0() << " Size: " << size << " GB\n";
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
            (find_permute->second)(permute.to_ta_expression());

        mpqc::time_point time1 = mpqc::now(world_, accurate_time_);
        double time = mpqc::duration_in_s(time0, time1);

        ExEnv::out0() << indent;
        ExEnv::out0() << "Permuted AO Integral: "
                      << utility::to_string(formula.string()) << " From "
                      << utility::to_string(permute.string());
        double size = mpqc::detail::array_size(result);
        ExEnv::out0() << " Size: " << size << " GB"
                      << " Time: " << time << " s\n";

        // store current array and delete old one
        ao_formula_registry_.insert(formula, result);
        //        ao_formula_registry_.purge_formula(world_,permute);
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

    result = ao_formula_registry_.retrieve(formula);
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
    ExEnv::out0() << indent;
    ExEnv::out0() << "Computed Inverse of Integral: "
                  << utility::to_string(formula.string());
    double size = mpqc::detail::array_size(result);

    time1 = mpqc::now(world_, accurate_time_);
    time += mpqc::duration_in_s(time0, time1);
    ExEnv::out0() << " Size: " << size << " GB"
                  << " Time: " << time << " s\n";
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

      ExEnv::out0() << indent;
      ExEnv::out0() << "Computed One Body Integral: "
                    << utility::to_string(formula.string());
      double size = mpqc::detail::array_size(result);
      ExEnv::out0() << " Size: " << size << " GB"
                    << " Time: " << time << " s\n";
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

      ExEnv::out0() << indent;
      ExEnv::out0() << "Computed Twobody Two Center Integral: "
                    << utility::to_string(formula.string());
      double size = mpqc::detail::array_size(result);
      ExEnv::out0() << " Size: " << size << " GB"
                    << " Time: " << time << " s\n";
    }
    // compute JK, requires orbital space registry
    else if (formula.oper().is_jk()) {
      // density fitting case

      // find the density
      auto space_index = detail::get_jk_orbital_space(formula.oper());
      auto& space = orbital_space_registry_->retrieve(space_index);

      auto obs = space.ao_index().name();
      if (formula.oper().has_option(Operator::Option::DensityFitting)) {
        auto three_center_formula = detail::get_jk_df_formula(formula, obs);

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
        auto space_index = detail::get_jk_orbital_space(formula.oper());
        auto& space = orbital_space_registry_->retrieve(space_index);
        auto obs = space.ao_index().name();
        // convert to ao formula
        auto four_center_formula = detail::get_jk_formula(formula, obs);
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

      ExEnv::out0() << indent;
      ExEnv::out0() << "Computed Coulumb/Exchange Integral: "
                    << utility::to_string(formula.string());
      double size = mpqc::detail::array_size(result);
      ExEnv::out0() << " Size: " << size << " GB"
                    << " Time: " << time << " s\n";
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

      ExEnv::out0() << indent;
      ExEnv::out0() << "Computed Coulumb/Exchange Integral: "
                    << utility::to_string(formula.string());
      double size = mpqc::detail::array_size(result);
      ExEnv::out0() << " Size: " << size << " GB"
                    << " Time: " << time << " s\n";
    }
    // compute Fock, requires orbital space registry
    else if (formula.oper().is_fock()) {
      auto formulas = detail::get_fock_formula(formula);

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

      ExEnv::out0() << indent;
      ExEnv::out0() << "Computed Fock Integral: "
                    << utility::to_string(formula.string());
      double size = mpqc::detail::array_size(result);
      ExEnv::out0() << " Size: " << size << " GB"
                    << " Time: " << time << " s\n";
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
    result = array_ops::eigen_to_array<Tile, Policy>(result.world(), result_eig,
                                                     tr_result, tr_result);

    if (formula.oper().type() == Operator::Type::cGTG ||
        formula.oper().type() == Operator::Type::cGTGCoulomb) {
      result("i,j") = -result("i,j");
    }
    time1 = mpqc::now(world_, accurate_time_);
    auto inv_time = mpqc::duration_in_s(time0, time1);
    ExEnv::out0() << indent;
    ExEnv::out0() << "Inverse Time: " << inv_time << " s\n";
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
      result = array_ops::eigen_to_array<Tile, Policy>(result.world(), inv_eig,
                                                       tr_result, tr_result);
    }

    if (formula.oper().type() == Operator::Type::cGTG ||
        formula.oper().type() == Operator::Type::cGTGCoulomb) {
      result("i,j") = -result("i,j");
    }
    time1 = mpqc::now(world_, accurate_time_);
    auto inv_sqrt_time = mpqc::duration_in_s(time0, time1);
    ExEnv::out0() << indent;
    ExEnv::out0() << "Inverse Square Root Time: " << inv_sqrt_time << " s\n";
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
  std::shared_ptr<Screener> p_screener = std::make_shared<Screener>(Screener{});

  parse_two_body_three_center(formula, engine_pool, bs_array, p_screener);
  result = compute_integrals(this->world_, engine_pool, bs_array, p_screener);

  time1 = mpqc::now(world_, accurate_time_);
  time += mpqc::duration_in_s(time0, time1);

  ExEnv::out0() << indent;
  ExEnv::out0() << "Computed Twobody Three Center Integral: "
                << utility::to_string(formula.string());
  double size = mpqc::detail::array_size(result);
  ExEnv::out0() << " Size: " << size << " GB"
                << " Time: " << time << " s\n";

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
    auto formula_strings = detail::get_df_formula(formula);

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

    ExEnv::out0() << indent;
    ExEnv::out0() << "Computed Twobody Four Center Density-Fitting Integral: "
                  << utility::to_string(formula.string());
    double size = mpqc::detail::array_size(result);
    ExEnv::out0() << " Size: " << size << " GB"
                  << " Time: " << time << " s\n";

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

    ExEnv::out0() << indent;
    ExEnv::out0() << "Computed Twobody Four Center Integral: "
                  << utility::to_string(formula.string());
    double size = mpqc::detail::array_size(result);
    ExEnv::out0() << " Size: " << size << " GB"
                  << " Time: " << time << " s\n";
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
  DirectTArray result;

  // find if in registry
  auto iter = direct_ao_formula_registry_.find(formula);

  if (iter != direct_ao_formula_registry_.end()) {
    result = iter->second;
    utility::print_par(world_, "Retrieved Direct AO Integral: ");
    utility::print_par(world_, utility::to_string(formula.string()));
    double size = mpqc::detail::array_size(result.array());
    utility::print_par(world_, " Size: ", size, " GB\n");
    return result;

  } else {
    // compute formula
    if (formula.rank() == 2) {
      result = compute_direct2(formula);
      direct_ao_formula_registry_.insert(formula, result);
    } else if (formula.rank() == 3) {
      result = compute_direct3(formula);
      direct_ao_formula_registry_.insert(formula, result);
    } else if (formula.rank() == 4) {
      result = compute_direct4(formula);
      direct_ao_formula_registry_.insert(formula, result);
    }
    return result;
  }
}

template <typename Tile, typename Policy>
typename AOFactory<Tile, Policy>::DirectTArray
AOFactory<Tile, Policy>::compute_direct2(const Formula& formula) {
  BasisVector bs_array;
  double time = 0.0;
  std::shared_ptr<utility::TSPool<libint2::Engine>> engine_pool;
  mpqc::time_point time0;
  mpqc::time_point time1;

  DirectTArray result;

  if (formula.oper().is_onebody()) {
    time0 = mpqc::now(world_, accurate_time_);

    parse_one_body(formula, engine_pool, bs_array);
    result = compute_direct_integrals(world_, engine_pool, bs_array);

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
    result = compute_direct_integrals(world_, engine_pool, bs_array);

    time1 = mpqc::now(world_, accurate_time_);
    time += mpqc::duration_in_s(time0, time1);

    utility::print_par(world_, "Computed Direct Twobody Two Center Integral: ");
    utility::print_par(world_, utility::to_string(formula.string()));
    double size = mpqc::detail::array_size(result.array());
    utility::print_par(world_, " Size: ", size, " GB");
    utility::print_par(world_, " Time: ", time, " s\n");
  } else {
    throw std::runtime_error("Unsupported Operator in DirectAOFactory!!\n");
  }

  madness::print_meminfo(
      world_.rank(), utility::wconcat("DirectAOFactory:", formula.string()));
  return result;
}

template <typename Tile, typename Policy>
typename AOFactory<Tile, Policy>::DirectTArray
AOFactory<Tile, Policy>::compute_direct3(const Formula& formula) {
  double time = 0.0;
  mpqc::time_point time0;
  mpqc::time_point time1;
  time0 = mpqc::now(world_, accurate_time_);
  DirectTArray result;

  BasisVector bs_array;
  std::shared_ptr<utility::TSPool<libint2::Engine>> engine_pool;
  std::shared_ptr<Screener> p_screener = std::make_shared<Screener>(Screener{});

  parse_two_body_three_center(formula, engine_pool, bs_array, p_screener);
  result =
      compute_direct_integrals(this->world_, engine_pool, bs_array, p_screener);

  time1 = mpqc::now(world_, accurate_time_);
  time += mpqc::duration_in_s(time0, time1);

  utility::print_par(world_, "Computed Direct Twobody Three Center Integral: ");
  utility::print_par(world_, utility::to_string(formula.string()));
  double size = mpqc::detail::array_size(result.array());
  utility::print_par(world_, " Size: ", size, " GB");
  utility::print_par(world_, " Time: ", time, " s\n");
  madness::print_meminfo(
      world_.rank(), utility::wconcat("DirectAOFactory:", formula.string()));

  return result;
}

template <typename Tile, typename Policy>
typename AOFactory<Tile, Policy>::DirectTArray
AOFactory<Tile, Policy>::compute_direct4(const Formula& formula) {
  if (formula.notation() != Formula::Notation::Chemical) {
    throw std::runtime_error(
        "Direct AO Integral Only Support Chemical Notation! \n");
  }

  double time = 0.0;
  mpqc::time_point time0;
  mpqc::time_point time1;
  DirectTArray result;

  time0 = mpqc::now(world_, accurate_time_);

  BasisVector bs_array;
  std::shared_ptr<Screener> p_screener = std::make_shared<Screener>(Screener{});
  std::shared_ptr<utility::TSPool<libint2::Engine>> engine_pool;

  parse_two_body_four_center(formula, engine_pool, bs_array, p_screener);

  result =
      compute_direct_integrals(this->world_, engine_pool, bs_array, p_screener);

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

  auto bra_basis =
      detail::index_to_basis(*this->orbital_basis_registry_, bra_index);
  auto ket_basis =
      detail::index_to_basis(*this->orbital_basis_registry_, ket_index);

  TA_ASSERT(bra_basis != nullptr);
  TA_ASSERT(ket_basis != nullptr);

  bases = BasisVector{{*bra_basis, *ket_basis}};

  auto oper_type = formula.oper().type();
  engine_pool = make_engine_pool(
      detail::to_libint2_operator(oper_type),
      utility::make_array_of_refs(*bra_basis, *ket_basis), libint2::BraKet::x_x,
      detail::to_libint2_operator_params(oper_type, *this->molecule_,
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

  auto bra_basis0 =
      detail::index_to_basis(*this->orbital_basis_registry_, bra_index0);
  auto ket_basis0 =
      detail::index_to_basis(*this->orbital_basis_registry_, ket_index0);

  TA_ASSERT(bra_basis0 != nullptr);
  TA_ASSERT(ket_basis0 != nullptr);

  bases = BasisVector{{*bra_basis0, *ket_basis0}};

  auto oper_type = formula.oper().type();
  engine_pool =
      make_engine_pool(detail::to_libint2_operator(oper_type),
                       utility::make_array_of_refs(*bra_basis0, *ket_basis0),
                       libint2::BraKet::xs_xs,
                       detail::to_libint2_operator_params(
                           oper_type, *this->molecule_, this->gtg_params_));
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

  auto bra_basis0 =
      detail::index_to_basis(*this->orbital_basis_registry_, bra_indexs[0]);
  auto ket_basis0 =
      detail::index_to_basis(*this->orbital_basis_registry_, ket_indexs[0]);
  auto ket_basis1 =
      detail::index_to_basis(*this->orbital_basis_registry_, ket_indexs[1]);

  TA_ASSERT(bra_basis0 != nullptr);
  TA_ASSERT(ket_basis0 != nullptr);
  TA_ASSERT(ket_basis1 != nullptr);

  bases = BasisVector{{*bra_basis0, *ket_basis0, *ket_basis1}};

  auto oper_type = formula.oper().type();
  engine_pool = make_engine_pool(
      detail::to_libint2_operator(oper_type),
      utility::make_array_of_refs(*bra_basis0, *ket_basis0, *ket_basis1),
      libint2::BraKet::xs_xx,
      detail::to_libint2_operator_params(oper_type, *this->molecule_,
                                         this->gtg_params_));

  if (!screen_.empty()) {
    /// make another engine to screener!!!

    auto screen_engine_pool = make_engine_pool(
        detail::to_libint2_operator(oper_type),
        utility::make_array_of_refs(*bra_basis0, *ket_basis0, *ket_basis1),
        libint2::BraKet::xx_xx,
        detail::to_libint2_operator_params(oper_type, *this->molecule_,
                                           this->gtg_params_));

    p_screener = detail::make_screener(this->world_, screen_engine_pool, bases,
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

  auto bra_basis0 =
      detail::index_to_basis(*this->orbital_basis_registry_, bra_indexs[0]);
  auto ket_basis0 =
      detail::index_to_basis(*this->orbital_basis_registry_, ket_indexs[0]);
  auto bra_basis1 =
      detail::index_to_basis(*this->orbital_basis_registry_, bra_indexs[1]);
  auto ket_basis1 =
      detail::index_to_basis(*this->orbital_basis_registry_, ket_indexs[1]);

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
      detail::to_libint2_operator_params(oper_type, *this->molecule_,
                                         this->gtg_params_));

  p_screener = detail::make_screener(this->world_, engine_pool, bases,
                                     this->screen_, this->screen_threshold_);
}

}  // namespace gaussian

}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_AO_FACTORY_IMPL_H_
