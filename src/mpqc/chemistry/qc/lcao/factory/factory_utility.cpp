//
// Created by Chong Peng on 3/2/16.
//

#include "mpqc/chemistry/qc/lcao/factory/factory_utility.h"

namespace mpqc {
namespace lcao {
namespace gaussian {

namespace detail {

libint2::Operator to_libint2_operator(Operator::Type mpqc_oper) {
  TA_USER_ASSERT((Operator::Type::__first_1body_operator <= mpqc_oper &&
                  mpqc_oper <= Operator::Type::__last_1body_operator) ||
                     (Operator::Type::__first_2body_operator <= mpqc_oper &&
                      mpqc_oper <= Operator::Type::__last_2body_operator),
                 "invalid Operator::Type");

  libint2::Operator result;
  switch (mpqc_oper) {
    case Operator::Type::Overlap: {
      result = libint2::Operator::overlap;
    } break;
    case Operator::Type::Kinetic: {
      result = libint2::Operator::kinetic;
    } break;
    case Operator::Type::Nuclear: {
      result = libint2::Operator::nuclear;
    } break;
    case Operator::Type::Coulomb: {
      result = libint2::Operator::coulomb;
    } break;
    case Operator::Type::Cadf: {
      result = libint2::Operator::coulomb;
    } break;
    case Operator::Type::cGTG: {
      result = libint2::Operator::cgtg;
    } break;
    case Operator::Type::cGTG2: {
      result = libint2::Operator::cgtg;
    } break;
    case Operator::Type::cGTGCoulomb: {
      result = libint2::Operator::cgtg_x_coulomb;
    } break;
    case Operator::Type::DelcGTG2: {
      result = libint2::Operator::delcgtg2;
    } break;
    default: {
      TA_USER_ASSERT(
          false, "mpqc::Operator::Type not convertible to libint2::Operator");
    }
  }
  return result;
}

libint2::any to_libint2_operator_params(
    Operator::Type mpqc_oper, const Molecule &molecule,
    const std::vector<std::pair<double, double>> &gtg_params) {
  TA_USER_ASSERT((Operator::Type::__first_1body_operator <= mpqc_oper &&
                  mpqc_oper <= Operator::Type::__last_1body_operator) ||
                     (Operator::Type::__first_2body_operator <= mpqc_oper &&
                      mpqc_oper <= Operator::Type::__last_2body_operator),
                 "invalid Operator::Type");

  libint2::any result;
  switch (mpqc_oper) {
    case Operator::Type::Nuclear: {
      result = make_q(molecule);
    } break;
    case Operator::Type::cGTG:
    case Operator::Type::cGTGCoulomb:
    case Operator::Type::DelcGTG2: {
      TA_USER_ASSERT(not gtg_params.empty(),
                     "Gaussian-type geminal not initialized");
      result = gtg_params;
    } break;
    case Operator::Type::cGTG2: {
      TA_USER_ASSERT(not gtg_params.empty(),
                     "Gaussian-type geminal not initialized");
      const auto &cgtg_params = gtg_params;
      const auto ng = cgtg_params.size();
      typename std::decay<decltype(cgtg_params)>::type cgtg2_params;
      cgtg2_params.reserve(ng * (ng + 1) / 2);
      for (auto b = 0; b < ng; ++b) {
        for (auto k = 0; k <= b; ++k) {
          const auto gexp = cgtg_params[b].first + cgtg_params[k].first;
          const auto gcoeff = cgtg_params[b].second * cgtg_params[k].second *
                              (b == k ? 1 : 2);  // if a != b include ab and ba
          cgtg2_params.push_back(std::make_pair(gexp, gcoeff));
        }
      }
      result = cgtg2_params;
    } break;
    default:;  // nothing to do
  }
  return result;
}

std::shared_ptr<Screener> make_screener(madness::World &world,
                                        const ShrPool<libint2::Engine> &engine,
                                        const BasisVector &bases,
                                        const std::string &screen,
                                        double screen_threshold) {
  std::shared_ptr<Screener> p_screen;
  if (screen.empty()) {
    p_screen = std::make_shared<Screener>(Screener{});
  } else if (screen == "qqr") {
    throw InputError("QQR screening is currently not availible", __FILE__,
                     __LINE__, "screen");
  } else if (screen == "schwarz") {
    p_screen = std::make_shared<gaussian::SchwarzScreen>(
        gaussian::create_schwarz_screener(world, engine, bases,
                                          screen_threshold));
  } else {
    throw InputError("Wrong screening method", __FILE__, __LINE__, "screen");
  }
  return p_screen;
}

std::shared_ptr<Basis> index_to_basis(
    const OrbitalBasisRegistry &basis_registry, const OrbitalIndex &index) {
  auto iter = basis_registry.find(index);
  if (iter == basis_registry.cend()) {
    throw std::runtime_error("Basis Set Not Found!!");
  } else {
    return iter->second;
  }
}

std::array<std::wstring, 3> get_df_formula(const Formula &formula) {
  std::array<std::wstring, 3> result;

  // chemical notation
  if (formula.notation() == Formula::Notation::Chemical) {
    std::wstring left = L"( Κ |" + formula.oper().oper_string() + L"| " +
                        formula.bra_indices()[0].name() + L" " +
                        formula.bra_indices()[1].name() + L" )";
    std::wstring right = L"( Κ |" + formula.oper().oper_string() + L"| " +
                         formula.ket_indices()[0].name() + L" " +
                         formula.ket_indices()[1].name() + L" )";
    std::wstring center =
        L"( Κ |" + formula.oper().oper_string() + L"| Λ)[inv]";
    result[0] = left;
    result[1] = center;
    result[2] = right;
  }
  // physical notation
  else {
    std::wstring left = L"( Κ |" + formula.oper().oper_string() + L"| " +
                        formula.bra_indices()[0].name() + L" " +
                        formula.ket_indices()[0].name() + L" )";
    std::wstring right = L"( Κ |" + formula.oper().oper_string() + L"| " +
                         formula.bra_indices()[1].name() + L" " +
                         formula.ket_indices()[1].name() + L" )";
    std::wstring center =
        L"( Κ |" + formula.oper().oper_string() + L"| Λ)[inv]";
    result[0] = left;
    result[1] = center;
    result[2] = right;
  }

  return result;
}

Formula get_jk_formula(const Formula &formula, const std::wstring &obs) {
  Formula result;

  Operator oper(L"G");
  result.set_operator(oper);
  if (formula.notation() == Formula::Notation::Chemical) {
    result.set_notation(Formula::Notation::Chemical);
    if (formula.oper().type() == Operator::Type::J) {
      result.bra_indices().push_back(formula.bra_indices()[0]);
      result.bra_indices().push_back(formula.ket_indices()[0]);
      result.ket_indices().push_back(obs + L"4");
      result.ket_indices().push_back(obs + L"5");

    } else {
      result.bra_indices().push_back(formula.bra_indices()[0]);
      result.bra_indices().push_back(obs + L"4");
      result.ket_indices().push_back(formula.ket_indices()[0]);
      result.ket_indices().push_back(obs + L"5");
    }
  } else {
    result.set_notation(Formula::Notation::Physical);
    if (formula.oper().type() == Operator::Type::K) {
      result.bra_indices().push_back(formula.bra_indices()[0]);
      result.bra_indices().push_back(formula.ket_indices()[0]);
      result.ket_indices().push_back(obs + L"4");
      result.ket_indices().push_back(obs + L"5");

    } else {
      result.bra_indices().push_back(formula.bra_indices()[0]);
      result.bra_indices().push_back(obs + L"4");
      result.ket_indices().push_back(formula.ket_indices()[0]);
      result.ket_indices().push_back(obs + L"5");
    }
  }
  return result;
}

std::array<Formula, 3> get_jk_df_formula(const Formula &formula,
                                         const std::wstring &obs) {
  std::array<Formula, 3> result;

  if (formula.oper().type() == Operator::Type::J) {
    std::wstring left = L"( Κ |G| " + formula.bra_indices()[0].name() + L" " +
                        formula.ket_indices()[0].name() + L" )";
    std::wstring right = L"( Κ |G| " + obs + L"4 " + obs + L"5 )";

    result[0] = Formula(left);
    result[2] = Formula(right);
  } else {
    std::wstring left =
        L"( Κ |G| " + formula.bra_indices()[0].name() + L" " + obs + L"4 )";
    std::wstring right =
        L"( Κ |G| " + formula.ket_indices()[0].name() + L" " + obs + L"5 )";

    result[0] = Formula(left);
    result[2] = Formula(right);
  }
  std::wstring center = L"( Κ |G| Λ)[inv]";
  result[1] = Formula(center);

  return result;
}

OrbitalIndex get_jk_orbital_space(const Operator &operation) {
  if (operation.type() == Operator::Type::J ||
      operation.type() == Operator::Type::K) {
    return OrbitalIndex(L"m");
  } else if (operation.type() == Operator::Type::KAlpha) {
    return OrbitalIndex(L"m_α");
  } else if (operation.type() == Operator::Type::KBeta) {
    return OrbitalIndex(L"m_β");
  } else {
    assert(false);
    return OrbitalIndex{};
  }
}

std::array<Formula, 3> get_fock_formula(const Formula &formula) {
  std::array<Formula, 3> result;
  Formula h(formula);
  h.set_operator(Operator(L"H"));
  decltype(h) j(formula);
  j.set_operator_type(Operator::Type::J);
  decltype(h) k(formula);
  if (formula.oper().type() == Operator::Type::Fock) {
    k.set_operator_type(Operator::Type::K);
  } else if (formula.oper().type() == Operator::Type::FockAlpha) {
    k.set_operator_type(Operator::Type::KAlpha);
  } else if (formula.oper().type() == Operator::Type::FockBeta) {
    k.set_operator_type(Operator::Type::KBeta);
  }

  result[0] = h;
  result[1] = j;
  result[2] = k;
  return result;
}

}  // namespace detail
}  // namespace gaussian

namespace detail {

bool if_all_lcao(const Formula &formula) {
  auto left = formula.bra_indices();
  for (auto &index : left) {
    if(!index.is_lcao()){
      return false;
    }
  }

  auto right = formula.ket_indices();
  for (auto &index : right) {
    if(!index.is_lcao()){
      return false;
    }
  }
  return true;
}

bool if_all_ao(const Formula &formula) {
  auto left = formula.bra_indices();
  for (auto &index : left) {
    if(!index.is_ao()){
      return false;
    }
  }

  auto right = formula.ket_indices();
  for (auto &index : right) {
    if(!index.is_ao()){
      return false;
    }
  }
  return true;
}

}  // namespace detail

}  // namespace lcao
}  // namespace mpqc
