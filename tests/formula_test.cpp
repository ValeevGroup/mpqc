//
// Created by Chong Peng on 10/19/15.
//

#include <catch.hpp>
#include <mpqc/chemistry/qc/expression/formula.h>

using namespace mpqc;

TEST_CASE("Formula Expression", "[formula]"){

    SECTION("one body test case"){

        Formula overlap(L"<κ|λ>");
        REQUIRE(overlap.oper().type() == Operator::Type::Overlap);
        REQUIRE(overlap.bra_indices().size() == 1);
        REQUIRE(overlap.ket_indices().size() == 1);
        REQUIRE(overlap.ket_indices()[0].index() == OrbitalIndex::Type::obs);
        REQUIRE(overlap.bra_indices()[0].index() == OrbitalIndex::Type::obs);
        REQUIRE(overlap.oper().is_onebody() == true);
        REQUIRE(overlap.notation() == Formula::Notation::Physical);
        REQUIRE(overlap.to_ta_expression() == "kappa, lamda");

        Formula J(L"<κ|J|λ>");
        REQUIRE( J.oper().type() == Operator::Type::J);
        REQUIRE( J.oper().is_jk() == true);

        Formula F(L"<κ|F(α)|λ>");
        REQUIRE( F.oper().type() == Operator::Type::FockAlpha);
        REQUIRE( F.oper().is_fock() == true);

    }

    SECTION("two body test case"){
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
        REQUIRE(three_center.to_ta_expression() == "KAPPA, kappa, lamda1");
    }

    SECTION("option test case"){
        Formula couloumb(L"<p q1|R|a' A'1> [df]");
        REQUIRE(couloumb.oper().has_option(Operator::Option::DensityFitting) == true);

        Formula couloumb1(L"<p q1|R|a' A'1> [df, inv_sqr]");
        Formula couloumb2(L"<p q1|R|a' A'1> [inv_sqr, df]");
        REQUIRE(couloumb1.oper().option() == couloumb2.oper().option());
    }

    SECTION("equality test case"){
        Formula formula1(L"<p q1|R|a' A'1>");
        Formula formula1df(L"<p q1|R|a' A'1>[df]");
        Formula formula2(L"<r s12|R|d' B'12>");
        Formula formula3(L"(r s12|R|d' B'12)");
        REQUIRE(formula1 == formula2);
        REQUIRE(formula2 != formula1df);
        REQUIRE(formula1 != formula3);
    }

    SECTION("comparision test case"){
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


    SECTION("error handling"){
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
