//
// Created by Chong Peng on 10/19/15.
//

#include "catch.hpp"

#include "mpqc/chemistry/qc/expression/orbital_index.h"

TEST_CASE("Orbital Index", "[orbital_index]") {
  SECTION("single letter case") {
    mpqc::OrbitalIndex m(L"m");
    REQUIRE(m.index() == mpqc::OrbitalIndex::Type::occ);
    REQUIRE(m.is_ao() == false);
    REQUIRE(m.is_mo() == true);
    REQUIRE(m.is_mo_in_obs() == true);
    REQUIRE(m.is_mo_in_abs() == false);
    REQUIRE(m.is_mo_in_ribs() == false);

    mpqc::OrbitalIndex n(L"n");
    REQUIRE(n.index() == mpqc::OrbitalIndex::Type::occ);
    REQUIRE(n.to_ta_expression() == "n");

    mpqc::OrbitalIndex i(L"i");
    REQUIRE(i.index() == mpqc::OrbitalIndex::Type::corr_occ);

    mpqc::OrbitalIndex l(L"l");
    REQUIRE(l.index() == mpqc::OrbitalIndex::Type::corr_occ);

    mpqc::OrbitalIndex x(L"x");
    REQUIRE(x.index() == mpqc::OrbitalIndex::Type::active);

    mpqc::OrbitalIndex y(L"y");
    REQUIRE(y.index() == mpqc::OrbitalIndex::Type::active);

    mpqc::OrbitalIndex a(L"a");
    REQUIRE(a.index() == mpqc::OrbitalIndex::Type::virt);

    mpqc::OrbitalIndex d(L"d");
    REQUIRE(d.index() == mpqc::OrbitalIndex::Type::virt);

    mpqc::OrbitalIndex p(L"p");
    REQUIRE(p.index() == mpqc::OrbitalIndex::Type::any);

    mpqc::OrbitalIndex s(L"p");
    REQUIRE(s.index() == mpqc::OrbitalIndex::Type::any);

    mpqc::OrbitalIndex k1(L"κ");
    REQUIRE(k1.index() == mpqc::OrbitalIndex::Type::obs);
    REQUIRE(k1.is_mo() == false);
    REQUIRE(k1.is_ao() == true);
    REQUIRE(k1.to_ta_expression() == "kappa");

    mpqc::OrbitalIndex v1(L"ν");
    REQUIRE(v1.index() == mpqc::OrbitalIndex::Type::obs);
    REQUIRE(v1.to_ta_expression() == "nu");

    mpqc::OrbitalIndex K1(L"Κ");
    REQUIRE(K1.index() == mpqc::OrbitalIndex::Type::dfbs);
    REQUIRE(K1.to_ta_expression() == "KAPPA");

    mpqc::OrbitalIndex V1(L"Ν");
    REQUIRE(V1.index() == mpqc::OrbitalIndex::Type::dfbs);
    REQUIRE(V1.to_ta_expression() == "NU");

    mpqc::OrbitalIndex a1(L"α");
    REQUIRE(a1.index() == mpqc::OrbitalIndex::Type::abs);
    REQUIRE(a1.to_ta_expression() == "alpha");

    mpqc::OrbitalIndex d1(L"δ");
    REQUIRE(d1.index() == mpqc::OrbitalIndex::Type::abs);
    REQUIRE(d1.to_ta_expression() == "delta");

    mpqc::OrbitalIndex p1(L"ρ");
    REQUIRE(p1.index() == mpqc::OrbitalIndex::Type::ribs);
    REQUIRE(p1.to_ta_expression() == "rho");

    mpqc::OrbitalIndex v2(L"υ");
    REQUIRE(v2.index() == mpqc::OrbitalIndex::Type::ribs);
    REQUIRE(v2.to_ta_expression() == "upsilon");
  }

  SECTION("letter with number") {
    mpqc::OrbitalIndex m(L"m2");
    REQUIRE(m.index() == mpqc::OrbitalIndex::Type::occ);
    REQUIRE(m.to_ta_expression() == "m2");

    mpqc::OrbitalIndex n(L"n6");
    REQUIRE(n.index() == mpqc::OrbitalIndex::Type::occ);

    mpqc::OrbitalIndex x(L"x7");
    REQUIRE(x.index() == mpqc::OrbitalIndex::Type::active);

    mpqc::OrbitalIndex y(L"y9");
    REQUIRE(y.index() == mpqc::OrbitalIndex::Type::active);

    mpqc::OrbitalIndex i(L"i1");
    REQUIRE(i.index() == mpqc::OrbitalIndex::Type::corr_occ);

    mpqc::OrbitalIndex l(L"l2");
    REQUIRE(l.index() == mpqc::OrbitalIndex::Type::corr_occ);

    mpqc::OrbitalIndex a(L"a14");
    REQUIRE(a.index() == mpqc::OrbitalIndex::Type::virt);

    mpqc::OrbitalIndex d(L"d25");
    REQUIRE(d.index() == mpqc::OrbitalIndex::Type::virt);

    mpqc::OrbitalIndex p(L"p36");
    REQUIRE(p.index() == mpqc::OrbitalIndex::Type::any);
    REQUIRE(p.to_ta_expression() == "p36");

    mpqc::OrbitalIndex s(L"s7");
    REQUIRE(s.index() == mpqc::OrbitalIndex::Type::any);

    mpqc::OrbitalIndex k1(L"κ1");
    REQUIRE(k1.index() == mpqc::OrbitalIndex::Type::obs);

    mpqc::OrbitalIndex v1(L"ν52");
    REQUIRE(v1.index() == mpqc::OrbitalIndex::Type::obs);

    mpqc::OrbitalIndex a1(L"α14");
    REQUIRE(a1.index() == mpqc::OrbitalIndex::Type::abs);
    REQUIRE(a1.to_ta_expression() == "alpha14");

    mpqc::OrbitalIndex d1(L"δ5");
    REQUIRE(d1.index() == mpqc::OrbitalIndex::Type::abs);
    REQUIRE(d1.to_ta_expression() == "delta5");

    mpqc::OrbitalIndex K1(L"Κ1");
    REQUIRE(K1.index() == mpqc::OrbitalIndex::Type::dfbs);
    REQUIRE(K1.to_ta_expression() == "KAPPA1");

    mpqc::OrbitalIndex V1(L"Ν2");
    REQUIRE(V1.index() == mpqc::OrbitalIndex::Type::dfbs);

    mpqc::OrbitalIndex p1(L"ρ123");
    REQUIRE(p1.index() == mpqc::OrbitalIndex::Type::ribs);

    mpqc::OrbitalIndex v2(L"υ45");
    REQUIRE(v2.index() == mpqc::OrbitalIndex::Type::ribs);
  }

  SECTION("letter with prime") {
    mpqc::OrbitalIndex m(L"m'");
    REQUIRE(m.index() == mpqc::OrbitalIndex::Type::frozen_occ);
    REQUIRE(m.to_ta_expression() == "m'");

    mpqc::OrbitalIndex n(L"n'");
    REQUIRE(n.index() == mpqc::OrbitalIndex::Type::frozen_occ);

    mpqc::OrbitalIndex P(L"P'");
    REQUIRE(P.index() == mpqc::OrbitalIndex::Type::allany);

    mpqc::OrbitalIndex S(L"S'");
    REQUIRE(S.index() == mpqc::OrbitalIndex::Type::allany);

    mpqc::OrbitalIndex a(L"a'");
    REQUIRE(a.index() == mpqc::OrbitalIndex::Type::othervirt);
    REQUIRE(a.is_mo_in_ribs() == false);
    REQUIRE(a.is_mo_in_abs() == true);
    REQUIRE(a.is_mo_in_obs() == false);

    mpqc::OrbitalIndex d(L"d'");
    REQUIRE(d.index() == mpqc::OrbitalIndex::Type::othervirt);

    mpqc::OrbitalIndex A(L"A'");
    REQUIRE(A.index() == mpqc::OrbitalIndex::Type::allvirt);
    REQUIRE(A.is_mo_in_ribs() == true);
    REQUIRE(A.is_mo_in_abs() == false);
    REQUIRE(A.is_mo_in_obs() == false);

    mpqc::OrbitalIndex D(L"D'");
    REQUIRE(D.index() == mpqc::OrbitalIndex::Type::allvirt);
  }

  SECTION("letter with number and prime") {
    mpqc::OrbitalIndex m(L"m'3");
    REQUIRE(m.index() == mpqc::OrbitalIndex::Type::frozen_occ);
    REQUIRE(m.to_ta_expression() == "m'3");

    mpqc::OrbitalIndex n(L"n'4");
    REQUIRE(n.index() == mpqc::OrbitalIndex::Type::frozen_occ);

    mpqc::OrbitalIndex P(L"P'1");
    REQUIRE(P.index() == mpqc::OrbitalIndex::Type::allany);

    mpqc::OrbitalIndex S(L"S'2");
    REQUIRE(S.index() == mpqc::OrbitalIndex::Type::allany);

    mpqc::OrbitalIndex a(L"a'4");
    REQUIRE(a.index() == mpqc::OrbitalIndex::Type::othervirt);

    mpqc::OrbitalIndex d(L"d'5");
    REQUIRE(d.index() == mpqc::OrbitalIndex::Type::othervirt);

    mpqc::OrbitalIndex A(L"A'16");
    REQUIRE(A.index() == mpqc::OrbitalIndex::Type::allvirt);

    mpqc::OrbitalIndex D(L"D'178");
    REQUIRE(D.index() == mpqc::OrbitalIndex::Type::allvirt);
  }

  SECTION("spin case") {
    mpqc::OrbitalIndex a(L"a'4_α");
    REQUIRE(a.spin() == mpqc::OrbitalIndex::Spin::Alpha);

    mpqc::OrbitalIndex m(L"λ");
    REQUIRE(m.spin() == mpqc::OrbitalIndex::Spin::None);
    REQUIRE(m.index() == mpqc::OrbitalIndex::Type::obs);

    mpqc::OrbitalIndex p(L"p_α");
    REQUIRE(p.spin() == mpqc::OrbitalIndex::Spin::Alpha);
    REQUIRE(p.index() == mpqc::OrbitalIndex::Type::any);

    mpqc::OrbitalIndex d1(L"δ5_β");
    REQUIRE(d1.spin() == mpqc::OrbitalIndex::Spin::Beta);
  }

  SECTION("comparison") {
    // spin
    mpqc::OrbitalIndex p(L"p_α");
    mpqc::OrbitalIndex q1(L"q23_β");
    mpqc::OrbitalIndex q2(L"r1");
    REQUIRE(p > q1);
    REQUIRE(p > q2);

    // index
    mpqc::OrbitalIndex A(L"A'16_β");
    mpqc::OrbitalIndex m(L"m'3");
    mpqc::OrbitalIndex k1(L"κ1_α");
    REQUIRE(m < A);
    REQUIRE(k1 < A);
    REQUIRE(k1 < m);
  }

  SECTION("error handling") {
    // key not allowed
    REQUIRE_THROWS(mpqc::OrbitalIndex(L"e"));
    REQUIRE_THROWS(mpqc::OrbitalIndex(L"ε"));
    REQUIRE_THROWS(mpqc::OrbitalIndex(L"I"));
    REQUIRE_THROWS(mpqc::OrbitalIndex(L"I'"));
    REQUIRE_THROWS(mpqc::OrbitalIndex(L"I1"));

    // wrong format
    //        REQUIRE_THROWS(mpqc::OrbitalIndex(L"1a3"));
    //        REQUIRE_THROWS(mpqc::OrbitalIndex(L"a4'"));

    // wrong spin
    REQUIRE_THROWS(mpqc::OrbitalIndex d1(L"δ5_a"));
    //        REQUIRE_THROWS(mpqc::OrbitalIndex d1(L"i_α  "));
  }
}