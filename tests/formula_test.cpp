//
// Created by Chong Peng on 10/19/15.
//

#include "catch.hpp"
#include <mpqc/chemistry/qc/expression/formula.h>
#include <mpqc/chemistry/qc/expression/permutation.h>

using namespace mpqc;

TEST_CASE("Formula Expression", "[formula]") {
  SECTION("one body test case") {
    Formula overlap(L"<κ|λ>");
    REQUIRE(overlap.oper().type() == Operator::Type::Overlap);
    REQUIRE(overlap.bra_indices().size() == 1);
    REQUIRE(overlap.ket_indices().size() == 1);
    REQUIRE(overlap.ket_indices()[0].index() == OrbitalIndex::Type::obs);
    REQUIRE(overlap.bra_indices()[0].index() == OrbitalIndex::Type::obs);
    REQUIRE(overlap.oper().is_onebody() == true);
    REQUIRE(overlap.notation() == Formula::Notation::Physical);
    //        REQUIRE(overlap.to_ta_expression() == "kappa, lamda");

    Formula J(L"<κ|J|λ>");
    REQUIRE(J.oper().type() == Operator::Type::J);
    REQUIRE(J.oper().is_jk() == true);

    Formula F(L"<κ|F(α)|λ>");
    REQUIRE(F.oper().type() == Operator::Type::FockAlpha);
    REQUIRE(F.oper().is_fock() == true);
  }

  SECTION("two body test case") {
    Formula kinetic(L"(p q|G|r s)");
    REQUIRE(kinetic.oper().type() == Operator::Type::Coulomb);
    REQUIRE(kinetic.notation() == Formula::Notation::Chemical);
    REQUIRE(kinetic.to_ta_expression() == "p, q, r, s");

    Formula r2(L"<p q| R2 |r s>");
    REQUIRE(r2.oper().type() == Operator::Type::cGTG2);
    Formula couloumb(L"<p_α q1|R|a' A'1>");

    REQUIRE(couloumb.bra_indices()[0].index() == OrbitalIndex::Type::any);
    REQUIRE(couloumb.bra_indices()[0].spin() == OrbitalIndex::Spin::Alpha);
    REQUIRE(couloumb.bra_indices()[1].index() == OrbitalIndex::Type::any);
    REQUIRE(couloumb.ket_indices()[0].index() == OrbitalIndex::Type::othervirt);
    REQUIRE(couloumb.ket_indices()[1].index() == OrbitalIndex::Type::allvirt);
    REQUIRE(couloumb.oper().is_twobody() == true);
    REQUIRE(couloumb.to_ta_expression() == "p_alpha, q1, a', A'1");

    Formula two_center(L"( a'_β |R|A'1)");
    REQUIRE(two_center.bra_indices().size() == 1);
    REQUIRE(two_center.to_ta_expression() == "a'_beta, A'1");

    Formula three_center(L"( Κ|R|κ λ1)");
    REQUIRE(three_center.bra_indices().size() == 1);
    REQUIRE(three_center.ket_indices().size() == 2);
    REQUIRE(three_center.to_ta_expression() == "KAPPA, kappa, lambda1");
  }

  SECTION("option test case") {
    Formula couloumb(L"<p q1|R|a' A'1> [df]");
    REQUIRE(couloumb.oper().has_option(Operator::Option::DensityFitting) ==
            true);

    Formula couloumb1(L"<p q1|R|a' A'1> [df, inv_sqr]");
    Formula couloumb2(L"<p q1|R|a' A'1> [inv_sqr, df]");
    REQUIRE(couloumb1.oper().option() == couloumb2.oper().option());
  }

  SECTION("equality test case") {
    Formula formula1(L"<p q1|R|a' A'1>");
    Formula formula1df(L"<p q1|R|a' A'1>[df]");
    Formula formula2(L"<r s12|R|d' B'12>");
    Formula formula3(L"(r s12|R|d' B'12)");
    REQUIRE(formula1 == formula2);
    REQUIRE(formula2 != formula1df);
    REQUIRE(formula1 != formula3);

    // special case, rank 2 doesn't compare notation
    Formula formula4(L"<i|H|a>");
    Formula formula5(L"(i|H|a)");
    REQUIRE(formula4 == formula5);
  }

  SECTION("comparision test case") {
    Formula formula1(L"<p q1|R|a' A'1>");
    Formula formula1df(L"<p q1|R|a' A'1>[df]");
    REQUIRE(formula1 < formula1df);

    Formula formula2(L"<r s12|G|d' B'12>");
    REQUIRE(formula2 < formula1);

    Formula formula3(L"<a b12|G|A' B'12>");
    REQUIRE(formula3 < formula2);

    Formula formula4(L"<p b12|G|d' B'12>[df]");
    REQUIRE(formula3 < formula4);
  }

  SECTION("permutation test") {
    // rank 2
    Formula formula1(L"<i|V|j>");
    auto permu1 = permutations(formula1);
    REQUIRE(permu1.empty());

    formula1 = Formula(L"<i|V|a>");
    permu1 = permutations(formula1);
    REQUIRE(permu1.size() == 1);
    REQUIRE(permu1[0] == Formula(L"<a|V|i>"));

    // rank 3
    Formula formula2(L"(Κ|G|a b)");
    auto permu2 = permutations(formula2);
    REQUIRE(permu2.empty());

    formula2 = Formula(L"(Κ|G|i a)");
    permu2 = permutations(formula2);
    REQUIRE(permu2.size() == 1);
    REQUIRE(permu2[0] == Formula(L"(Κ|G|a i)"));

    // rank4
    Formula formula3(L"(i j|G|a b)");
    auto permu3 = permutations(formula3);
    REQUIRE(permu3.size() == 1);
    REQUIRE(permu3[0] == Formula(L"(a b|G|i j)"));

    formula3 = Formula(L"<i j|G|a b>");
    permu3 = permutations(formula3);
    REQUIRE(permu3.size() == 3);
    REQUIRE(std::find(permu3.cbegin(), permu3.cend(),
                      Formula(L"<a j|G|i b>")) != permu3.cend());
    REQUIRE(std::find(permu3.cbegin(), permu3.cend(),
                      Formula(L"<a b|G|i j>")) != permu3.cend());
    REQUIRE(std::find(permu3.cbegin(), permu3.cend(),
                      Formula(L"<i b|G|a j>")) != permu3.cend());

    formula3 = Formula(L"(i a|G|p b)");
    permu3 = permutations(formula3);
    REQUIRE(permu3.size() == 7);
    REQUIRE(std::find(permu3.cbegin(), permu3.cend(),
                      Formula(L"(a i|G|p b)")) != permu3.cend());
    REQUIRE(std::find(permu3.cbegin(), permu3.cend(),
                      Formula(L"(i a|G|b p)")) != permu3.cend());
    REQUIRE(std::find(permu3.cbegin(), permu3.cend(),
                      Formula(L"(a i|G|b p)")) != permu3.cend());
    REQUIRE(std::find(permu3.cbegin(), permu3.cend(),
                      Formula(L"(p b|G|i a)")) != permu3.cend());
    REQUIRE(std::find(permu3.cbegin(), permu3.cend(),
                      Formula(L"(b p|G|i a)")) != permu3.cend());
    REQUIRE(std::find(permu3.cbegin(), permu3.cend(),
                      Formula(L"(p b|G|a i)")) != permu3.cend());
    REQUIRE(std::find(permu3.cbegin(), permu3.cend(),
                      Formula(L"(b p|G|a i)")) != permu3.cend());

    formula3 = Formula(L"<i p|G|a b>");
    permu3 = permutations(formula3);
    REQUIRE(permu3.size() == 7);

    formula3 = Formula(L"<i p|G| a b>");
    permu3 = permutations(formula3);
    REQUIRE(permu3.size() == 7);
  }

  SECTION("error handling") {
    // wrong operation
    REQUIRE_THROWS(Formula kinetic(L"<p q|t|r s>"));
    //
    //        // wrong format
    REQUIRE_THROWS(Formula kinetic(L"<p q|T|r s>"));
    REQUIRE_THROWS(Formula kinetic(L"<p |q|K|r s>"));
    //
    REQUIRE_THROWS(Formula kinetic(L"<pq|T|r s>"));
    //        REQUIRE_THROWS(Formula kinetic(L"<p1q2|T|r s>"));
    //        REQUIRE_THROWS(Formula kinetic(L"<a'q2|T|r s>"));
  }
}
