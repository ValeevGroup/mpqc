//
// Created by Chong Peng on 10/19/15.
//

#include <catch.hpp>
#include "../expression/formula.h"

using namespace mpqc;

TEST_CASE("Formula Expression", "[formula]"){

    SECTION("operation test case"){

        Formula overlap(L"<a|b>");
        REQUIRE(overlap.operation() == Formula::Operation::Overlap);

        Formula kinetic(L"<p q|T|r s>");
        REQUIRE(kinetic.operation() == Formula::Operation::Kinetic);

        Formula r2(L"<p q| R2 |r s>");
        REQUIRE(r2.operation() == Formula::Operation::cGTG2);
    }

    SECTION("left right index case"){
        Formula kinetic(L"<p q1|T|a' A1'>");
        REQUIRE(kinetic.left_index()[0].index() == OrbitalIndex::Index::any);
        REQUIRE(kinetic.left_index()[1].index() == OrbitalIndex::Index::any);
        REQUIRE(kinetic.right_index()[0].index() == OrbitalIndex::Index::othervirt);
        REQUIRE(kinetic.right_index()[1].index() == OrbitalIndex::Index::allvirt);

    }

    SECTION("error handling"){
        // wrong operation
        REQUIRE_THROWS(Formula kinetic(L"<p q|t|r s>"));

        // wrong format
        REQUIRE_THROWS(Formula kinetic(L"p q|T|r s>"));
        REQUIRE_THROWS(Formula kinetic(L"<p |q|K|r s>"));

        REQUIRE_THROWS(Formula kinetic(L"<pq|T|r s>"));
        REQUIRE_THROWS(Formula kinetic(L"<p1q2|T|r s>"));
        REQUIRE_THROWS(Formula kinetic(L"<a'q2|T|r s>"));

    }

}
