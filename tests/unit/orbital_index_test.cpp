//
// Created by Chong Peng on 10/19/15.
//

#include "catch.hpp"

#include "mpqc/chemistry/qc/lcao/expression/orbital_index.h"

using mpqc::lcao::OrbitalIndex;

TEST_CASE("Orbital Index", "[orbital_index]") {
  SECTION("single letter case") {
    OrbitalIndex m(L"m");
    REQUIRE(m.index() == OrbitalIndex::Type::occ);
    REQUIRE(m.is_ao() == false);
    REQUIRE(m.is_mo() == true);
    REQUIRE(m.is_mo_in_obs() == true);
    REQUIRE(m.is_mo_in_abs() == false);
    REQUIRE(m.is_mo_in_ribs() == false);

    OrbitalIndex n(L"n");
    REQUIRE(n.index() == OrbitalIndex::Type::occ);
    REQUIRE(n.to_ta_expression() == "n");

    OrbitalIndex i(L"i");
    REQUIRE(i.index() == OrbitalIndex::Type::corr_occ);

    OrbitalIndex l(L"l");
    REQUIRE(l.index() == OrbitalIndex::Type::corr_occ);

    OrbitalIndex x(L"x");
    REQUIRE(x.index() == OrbitalIndex::Type::active);

    OrbitalIndex y(L"y");
    REQUIRE(y.index() == OrbitalIndex::Type::active);

    OrbitalIndex a(L"a");
    REQUIRE(a.index() == OrbitalIndex::Type::virt);

    OrbitalIndex d(L"d");
    REQUIRE(d.index() == OrbitalIndex::Type::virt);

    OrbitalIndex p(L"p");
    REQUIRE(p.index() == OrbitalIndex::Type::any);

    OrbitalIndex s(L"p");
    REQUIRE(s.index() == OrbitalIndex::Type::any);

    OrbitalIndex k1(L"κ");
    REQUIRE(k1.index() == OrbitalIndex::Type::obs);
    REQUIRE(k1.is_mo() == false);
    REQUIRE(k1.is_ao() == true);
    REQUIRE(k1.to_ta_expression() == "kappa");

    OrbitalIndex v1(L"ν");
    REQUIRE(v1.index() == OrbitalIndex::Type::obs);
    REQUIRE(v1.to_ta_expression() == "nu");

    OrbitalIndex K1(L"Κ");
    REQUIRE(K1.index() == OrbitalIndex::Type::dfbs);
    REQUIRE(K1.to_ta_expression() == "KAPPA");

    OrbitalIndex V1(L"Ν");
    REQUIRE(V1.index() == OrbitalIndex::Type::dfbs);
    REQUIRE(V1.to_ta_expression() == "NU");

    OrbitalIndex a1(L"α");
    REQUIRE(a1.index() == OrbitalIndex::Type::abs);
    REQUIRE(a1.to_ta_expression() == "alpha");

    OrbitalIndex d1(L"δ");
    REQUIRE(d1.index() == OrbitalIndex::Type::abs);
    REQUIRE(d1.to_ta_expression() == "delta");

    OrbitalIndex p1(L"ρ");
    REQUIRE(p1.index() == OrbitalIndex::Type::ribs);
    REQUIRE(p1.to_ta_expression() == "rho");

    OrbitalIndex v2(L"υ");
    REQUIRE(v2.index() == OrbitalIndex::Type::ribs);
    REQUIRE(v2.to_ta_expression() == "upsilon");
  }

  SECTION("letter with number") {
    OrbitalIndex m(L"m2");
    REQUIRE(m.index() == OrbitalIndex::Type::occ);
    REQUIRE(m.to_ta_expression() == "m2");

    OrbitalIndex n(L"n6");
    REQUIRE(n.index() == OrbitalIndex::Type::occ);

    OrbitalIndex x(L"x7");
    REQUIRE(x.index() == OrbitalIndex::Type::active);

    OrbitalIndex y(L"y9");
    REQUIRE(y.index() == OrbitalIndex::Type::active);

    OrbitalIndex i(L"i1");
    REQUIRE(i.index() == OrbitalIndex::Type::corr_occ);

    OrbitalIndex l(L"l2");
    REQUIRE(l.index() == OrbitalIndex::Type::corr_occ);

    OrbitalIndex a(L"a14");
    REQUIRE(a.index() == OrbitalIndex::Type::virt);

    OrbitalIndex d(L"d25");
    REQUIRE(d.index() == OrbitalIndex::Type::virt);

    OrbitalIndex p(L"p36");
    REQUIRE(p.index() == OrbitalIndex::Type::any);
    REQUIRE(p.to_ta_expression() == "p36");

    OrbitalIndex s(L"s7");
    REQUIRE(s.index() == OrbitalIndex::Type::any);

    OrbitalIndex k1(L"κ1");
    REQUIRE(k1.index() == OrbitalIndex::Type::obs);

    OrbitalIndex v1(L"ν52");
    REQUIRE(v1.index() == OrbitalIndex::Type::obs);

    OrbitalIndex a1(L"α14");
    REQUIRE(a1.index() == OrbitalIndex::Type::abs);
    REQUIRE(a1.to_ta_expression() == "alpha14");

    OrbitalIndex d1(L"δ5");
    REQUIRE(d1.index() == OrbitalIndex::Type::abs);
    REQUIRE(d1.to_ta_expression() == "delta5");

    OrbitalIndex K1(L"Κ1");
    REQUIRE(K1.index() == OrbitalIndex::Type::dfbs);
    REQUIRE(K1.to_ta_expression() == "KAPPA1");

    OrbitalIndex V1(L"Ν2");
    REQUIRE(V1.index() == OrbitalIndex::Type::dfbs);

    OrbitalIndex p1(L"ρ123");
    REQUIRE(p1.index() == OrbitalIndex::Type::ribs);

    OrbitalIndex v2(L"υ45");
    REQUIRE(v2.index() == OrbitalIndex::Type::ribs);
  }

  SECTION("letter with prime") {
    OrbitalIndex m(L"m'");
    REQUIRE(m.index() == OrbitalIndex::Type::frozen_occ);
    REQUIRE(m.to_ta_expression() == "m'");

    OrbitalIndex n(L"n'");
    REQUIRE(n.index() == OrbitalIndex::Type::frozen_occ);

    OrbitalIndex P(L"P'");
    REQUIRE(P.index() == OrbitalIndex::Type::allany);

    OrbitalIndex S(L"S'");
    REQUIRE(S.index() == OrbitalIndex::Type::allany);

    OrbitalIndex a(L"a'");
    REQUIRE(a.index() == OrbitalIndex::Type::othervirt);
    REQUIRE(a.is_mo_in_ribs() == false);
    REQUIRE(a.is_mo_in_abs() == true);
    REQUIRE(a.is_mo_in_obs() == false);

    OrbitalIndex d(L"d'");
    REQUIRE(d.index() == OrbitalIndex::Type::othervirt);

    OrbitalIndex A(L"A'");
    REQUIRE(A.index() == OrbitalIndex::Type::allvirt);
    REQUIRE(A.is_mo_in_ribs() == true);
    REQUIRE(A.is_mo_in_abs() == false);
    REQUIRE(A.is_mo_in_obs() == false);

    OrbitalIndex D(L"D'");
    REQUIRE(D.index() == OrbitalIndex::Type::allvirt);
  }

  SECTION("letter with number and prime") {
    OrbitalIndex m(L"m'3");
    REQUIRE(m.index() == OrbitalIndex::Type::frozen_occ);
    REQUIRE(m.to_ta_expression() == "m'3");

    OrbitalIndex n(L"n'4");
    REQUIRE(n.index() == OrbitalIndex::Type::frozen_occ);

    OrbitalIndex P(L"P'1");
    REQUIRE(P.index() == OrbitalIndex::Type::allany);

    OrbitalIndex S(L"S'2");
    REQUIRE(S.index() == OrbitalIndex::Type::allany);

    OrbitalIndex a(L"a'4");
    REQUIRE(a.index() == OrbitalIndex::Type::othervirt);

    OrbitalIndex d(L"d'5");
    REQUIRE(d.index() == OrbitalIndex::Type::othervirt);

    OrbitalIndex A(L"A'16");
    REQUIRE(A.index() == OrbitalIndex::Type::allvirt);

    OrbitalIndex D(L"D'178");
    REQUIRE(D.index() == OrbitalIndex::Type::allvirt);
  }

  SECTION("spin case") {
    OrbitalIndex a(L"a'4_α");
    REQUIRE(a.spin() == OrbitalIndex::Spin::Alpha);

    OrbitalIndex m(L"λ");
    REQUIRE(m.spin() == OrbitalIndex::Spin::None);
    REQUIRE(m.index() == OrbitalIndex::Type::obs);

    OrbitalIndex p(L"p_α");
    REQUIRE(p.spin() == OrbitalIndex::Spin::Alpha);
    REQUIRE(p.index() == OrbitalIndex::Type::any);

    OrbitalIndex d1(L"δ5_β");
    REQUIRE(d1.spin() == OrbitalIndex::Spin::Beta);
  }

  SECTION("comparison") {
    // spin
    OrbitalIndex p(L"p_α");
    OrbitalIndex q1(L"q23_β");
    OrbitalIndex q2(L"r1");
    REQUIRE(p > q1);
    REQUIRE(p > q2);

    // index
    OrbitalIndex A(L"A'16_β");
    OrbitalIndex m(L"m'3");
    OrbitalIndex k1(L"κ1_α");
    REQUIRE(m < A);
    REQUIRE(k1 < A);
    REQUIRE(k1 < m);
  }

  SECTION("error handling") {
    // key not allowed
    REQUIRE_THROWS(OrbitalIndex(L"e"));
    REQUIRE_THROWS(OrbitalIndex(L"ε"));
    REQUIRE_THROWS(OrbitalIndex(L"I"));
    REQUIRE_THROWS(OrbitalIndex(L"I'"));
    REQUIRE_THROWS(OrbitalIndex(L"I1"));

    // wrong format
    //        REQUIRE_THROWS(OrbitalIndex(L"1a3"));
    //        REQUIRE_THROWS(OrbitalIndex(L"a4'"));

    // wrong spin
    REQUIRE_THROWS(OrbitalIndex d1(L"δ5_a"));
    //        REQUIRE_THROWS(OrbitalIndex d1(L"i_α  "));
  }
}
