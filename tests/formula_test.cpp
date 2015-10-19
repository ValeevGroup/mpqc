//
// Created by Chong Peng on 10/19/15.
//

#include <catch.hpp>
#include "../expression/formula.h"

using namespace mpqc;

TEST_CASE("Formula Expression", "[formula]"){

    SECTION("operation test case"){
        Formula kinetic("<p q|K|r s>");
        REQUIRE(kinetic.operation() == Formula::Operation::Kinetic);

        Formula r2("<p q|R2|r s>");
        REQUIRE(r2.operation() == Formula::Operation::cGTG2);
    }

    SECTION("left right index case"){
        Formula kinetic("<p q1|K|a' A1'>");
        REQUIRE(kinetic.left_index()[0].index() == OrbitalIndex::Index::any);
        REQUIRE(kinetic.left_index()[1].index() == OrbitalIndex::Index::any);
        REQUIRE(kinetic.left_index()[2].index() == OrbitalIndex::Index::othervirt);
        REQUIRE(kinetic.left_index()[3].index() == OrbitalIndex::Index::allvirt);

    }

    SECTION("error handling"){
        // wrong operation
        REQUIRE_THROWS(Formula kinetic("<p q|k|r s>"));

        // wrong format
        REQUIRE_THROWS(Formula kinetic("p q|K|r s>"));

        REQUIRE_THROWS(Formula kinetic("<pq|K|r s>"));
        REQUIRE_THROWS(Formula kinetic("<p1q2|K|r s>"));
        REQUIRE_THROWS(Formula kinetic("<a'q2|K|r s>"));

    }

}
