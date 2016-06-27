//
// Created by Chong Peng on 3/2/16.
//

#include <mpqc/chemistry/qc/integrals/atomic_integral_base.h>

namespace mpqc{
namespace integrals{

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

libint2::any to_libint2_operator_params(Operator::Type mpqc_oper,
                                        const AtomicIntegralBase& base) {
  TA_USER_ASSERT((Operator::Type::__first_1body_operator <= mpqc_oper &&
                  mpqc_oper <= Operator::Type::__last_1body_operator) ||
                 (Operator::Type::__first_2body_operator <= mpqc_oper &&
                  mpqc_oper <= Operator::Type::__last_2body_operator),
                  "invalid Operator::Type");

  libint2::any result;
  switch (mpqc_oper) {
    case Operator::Type::Nuclear: {
      result = make_q(base.molecule());
    } break;
    case Operator::Type::cGTG:
    case Operator::Type::cGTGCoulomb:
    case Operator::Type::DelcGTG2:
    {
      result = base.gtg_params();
    } break;
    case Operator::Type::cGTG2: {
      const auto& cgtg_params = base.gtg_params();
      const auto ng = cgtg_params.size();
      std::decay<decltype(cgtg_params)>::type cgtg2_params;
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
    default: ;  // nothing to do
  }
  return result;
}

libint2::Engine AtomicIntegralBase::make_engine(const Operator& oper,
                                                int64_t max_nprim,
                                                int64_t max_am) {
  auto op = to_libint2_operator(oper.type());
  auto params = to_libint2_operator_params(oper.type(), *this);
  libint2::Engine engine(op, max_nprim, static_cast<int>(max_am), 0);
  engine.set_params(std::move(params));

  return engine;
}

void AtomicIntegralBase::parse_one_body(const Formula& formula,
                                        std::shared_ptr<EnginePool<libint2::Engine>>& engine_pool,
                                        Barray<2>& bases){

    auto bra_indices = formula.left_index();
    auto ket_indices = formula.right_index();

    TA_ASSERT(bra_indices.size() == 1);
    TA_ASSERT(ket_indices.size() == 1);

    auto bra_index = bra_indices[0];
    auto ket_index = ket_indices[0];

    TA_ASSERT(bra_index.is_ao());
    TA_ASSERT(ket_index.is_ao());

    auto bra_basis = this->index_to_basis(bra_index);
    auto ket_basis = this->index_to_basis(ket_index);

    TA_ASSERT(bra_basis != nullptr);
    TA_ASSERT(ket_basis != nullptr);

    auto max_nprim = std::max(bra_basis->max_nprim(), ket_basis->max_nprim());
    auto max_am = std::max(bra_basis->max_am(), ket_basis->max_am());
    bases = mpqc::utility::make_array(*bra_basis, *ket_basis);

    engine_pool = make_pool(make_engine(formula.oper(),max_nprim,max_am));
}

void AtomicIntegralBase::parse_two_body_two_center(
    const Formula &formula,
    std::shared_ptr<EnginePool<libint2::Engine>> &engine_pool,
    Barray<2> &bases) {

    TA_USER_ASSERT(formula.notation() == Formula::Notation::Chemical, "Two Body Two Center Integral Must Use Chemical Notation");

    auto bra_indexs = formula.left_index();
    auto ket_indexs = formula.right_index();

    TA_ASSERT(bra_indexs.size() == 1);
    TA_ASSERT(ket_indexs.size() == 1);

    auto bra_index0 = bra_indexs[0];
    auto ket_index0 = ket_indexs[0];

    TA_ASSERT(ket_index0.is_ao());
    TA_ASSERT(bra_index0.is_ao());

    auto bra_basis0 = this->index_to_basis(bra_index0);
    auto ket_basis0 = this->index_to_basis(ket_index0);

    TA_ASSERT(bra_basis0 != nullptr);
    TA_ASSERT(ket_basis0 != nullptr);

    auto max_nprim = std::max(bra_basis0->max_nprim(), ket_basis0->max_nprim());
    auto max_am = std::max(bra_basis0->max_am(), ket_basis0->max_am());
    bases = mpqc::utility::make_array(*bra_basis0, *ket_basis0);

    engine_pool = make_pool(make_engine(formula.oper(),max_nprim,max_am));
}

void AtomicIntegralBase::parse_two_body_three_center(
    const Formula &formula,
    std::shared_ptr<EnginePool<libint2::Engine>> &engine_pool,
    Barray<3> &bases) {

    TA_USER_ASSERT(formula.notation() == Formula::Notation::Chemical, "Three Center Integral Must Use Chemical Notation");

    auto bra_indexs = formula.left_index();
    auto ket_indexs = formula.right_index();

    TA_ASSERT(bra_indexs.size() == 1);
    TA_ASSERT(ket_indexs.size() == 2);


    TA_ASSERT(bra_indexs[0].is_ao());
    TA_ASSERT(ket_indexs[0].is_ao());
    TA_ASSERT(ket_indexs[1].is_ao());

    auto bra_basis0 = this->index_to_basis(bra_indexs[0]);
    auto ket_basis0 = this->index_to_basis(ket_indexs[0]);
    auto ket_basis1 = this->index_to_basis(ket_indexs[1]);

    TA_ASSERT(bra_basis0 != nullptr);
    TA_ASSERT(ket_basis0 != nullptr);
    TA_ASSERT(ket_basis1 != nullptr);

    auto max_nprim = std::max({bra_basis0->max_nprim(), ket_basis0->max_nprim()
                                 ,ket_basis1->max_nprim()});
    auto max_am = std::max({bra_basis0->max_am(), ket_basis0->max_am(),
                       ket_basis1->max_am()});

    bases = mpqc::utility::make_array(*bra_basis0, *ket_basis0, *ket_basis1);

    engine_pool = make_pool(make_engine(formula.oper(),max_nprim,max_am));
}

void AtomicIntegralBase::parse_two_body_four_center(
    const Formula &formula,
    std::shared_ptr<EnginePool<libint2::Engine>> &engine_pool,
    Barray<4> &bases) {

    auto bra_indexs = formula.left_index();
    auto ket_indexs = formula.right_index();

    TA_ASSERT(bra_indexs.size() == 2);
    TA_ASSERT(ket_indexs.size() == 2);


    TA_ASSERT(bra_indexs[0].is_ao());
    TA_ASSERT(ket_indexs[0].is_ao());
    TA_ASSERT(bra_indexs[1].is_ao());
    TA_ASSERT(ket_indexs[1].is_ao());

    auto bra_basis0 = this->index_to_basis(bra_indexs[0]);
    auto ket_basis0 = this->index_to_basis(ket_indexs[0]);
    auto bra_basis1 = this->index_to_basis(bra_indexs[1]);
    auto ket_basis1 = this->index_to_basis(ket_indexs[1]);

    TA_ASSERT(bra_basis0 != nullptr);
    TA_ASSERT(ket_basis0 != nullptr);
    TA_ASSERT(bra_basis1 != nullptr);
    TA_ASSERT(ket_basis1 != nullptr);

    auto max_nprim = std::max({bra_basis0->max_nprim(), bra_basis1->max_nprim()
                                 ,ket_basis0->max_nprim(), ket_basis1->max_nprim()});
    auto max_am = std::max({bra_basis0->max_am(), bra_basis1->max_am(),
                       ket_basis0->max_am(),ket_basis1->max_am()});

    if (formula.notation() == Formula::Notation::Chemical){
        bases = {{*bra_basis0, *bra_basis1, *ket_basis0, *ket_basis1}};
    }
    else{
        bases = {{*bra_basis0, *ket_basis0, *bra_basis1, *ket_basis1}};
    }

    engine_pool = make_pool(make_engine(formula.oper(),max_nprim,max_am));
}


std::array<std::wstring, 3> AtomicIntegralBase::get_df_formula(const Formula &formula) {
    std::array<std::wstring, 3> result;

    //chemical notation
    if (formula.notation() == Formula::Notation::Chemical) {

        std::wstring left =
                L"( Κ |" + formula.oper().oper_string() + L"| " + formula.left_index()[0].name() + L" " +
                formula.left_index()[1].name() + L" )";
        std::wstring right =
                L"( Κ |" + formula.oper().oper_string() + L"| " + formula.right_index()[0].name() + L" " +
                formula.right_index()[1].name() + L" )";
        std::wstring center = L"( Κ |" + formula.oper().oper_string() + L"| Λ)[inv]";
        result[0] = left;
        result[1] = center;
        result[2] = right;
    }
    //physical notation
    else {
        std::wstring left =
                L"( Κ |" + formula.oper().oper_string() + L"| " + formula.left_index()[0].name() + L" " +
                formula.right_index()[0].name() + L" )";
        std::wstring right =
                L"( Κ |" + formula.oper().oper_string() + L"| " + formula.left_index()[1].name() + L" " +
                formula.right_index()[1].name() + L" )";
        std::wstring center = L"( Κ |" + formula.oper().oper_string() + L"| Λ)[inv]";
        result[0] = left;
        result[1] = center;
        result[2] = right;
    }

    return result;
}


Formula AtomicIntegralBase::get_jk_formula(const Formula &formula, const std::wstring& obs) {

    Formula result;

    result.set_notation(Formula::Notation::Chemical);
    Operator oper(L"G");
    result.set_operator(oper);
    if(formula.oper().type() == Operator::Type::J){

        result.left_index().push_back(formula.left_index()[0]);
        result.left_index().push_back(formula.right_index()[0]);
        result.right_index().push_back(obs);
        result.right_index().push_back(obs);

    }
    else{

        result.left_index().push_back(formula.left_index()[0]);
        result.left_index().push_back(obs);
        result.right_index().push_back(formula.right_index()[0]);
        result.right_index().push_back(obs);

    }
    return result;
}

std::array<Formula,3> AtomicIntegralBase::get_jk_df_formula(const Formula &formula, const std::wstring& obs) {
    std::array<Formula,3> result;

    if(formula.oper().type() == Operator::Type::J){
        std::wstring left =  L"( Κ |G| " + formula.left_index()[0].name() + L" " + formula.right_index()[0].name() + L" )";
        std::wstring right = L"( Κ |G| " + obs + L" " + obs + L" )";

        result[0] = Formula(left);
        result[2] = Formula(right);
    }
    else{
        std::wstring left =  L"( Κ |G| " + formula.left_index()[0].name() + L" " + obs  + L" )";
        std::wstring right = L"( Κ |G| " + formula.right_index()[0].name() +L" " + obs + L" )";

        result[0] = Formula(left);
        result[2] = Formula(right);
    }
    std::wstring center = L"( Κ |G| Λ)[inv]";
    result[1] = Formula(center);

    return result;
}

OrbitalIndex AtomicIntegralBase::get_jk_orbital_space(const Operator &operation) {

    if(operation.type() == Operator::Type::J || operation.type() == Operator::Type::K){
        return OrbitalIndex(L"m");
    }
    else if(operation.type() == Operator::Type::KAlpha){
        return OrbitalIndex(L"m_α");
    }
    else if(operation.type() == Operator::Type::KBeta){
        return OrbitalIndex(L"m_β");
    } else {
        assert(false);
        return OrbitalIndex{};
    }
}

std::array<Formula, 3> AtomicIntegralBase::get_fock_formula(const Formula &formula) {

    std::array<Formula, 3> result;
    Formula h(formula);
    h.set_operator(Operator(L"H"));
    decltype(h) j(formula);
    j.set_operator_type(Operator::Type::J);
    decltype(h) k(formula);
    if(formula.oper().type() == Operator::Type::Fock){
      k.set_operator_type(Operator::Type::K);
    }
    else if(formula.oper().type() == Operator::Type::FockAlpha){
      k.set_operator_type(Operator::Type::KAlpha);
    }
    else if(formula.oper().type() == Operator::Type::FockBeta){
      k.set_operator_type(Operator::Type::KBeta);
    }

    result[0] = h;
    result[1] = j;
    result[2] = k;
    return result;
}
}// end of namespace integral
}// end of namespace mpqc
