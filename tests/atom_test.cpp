#include "catch.hpp"

#include "../molecule/atom.h"

SCENARIO("atoms can be intialized", "[atom]"){
    GIVEN("a default initalized atom"){
        tcc::molecule::Atom a;
        auto center = a.center();
        REQUIRE(center[0] == 0);
        REQUIRE(center[1] == 0);
        REQUIRE(center[2] == 0);

        REQUIRE(a.mass() == 0);
        REQUIRE(a.charge() == 0);

        WHEN("the atom is assigned to a new atom"){
            a = tcc::molecule::Atom({1.0, 2.0, 3.0}, 1.07, 1);
            center = a.center();
            REQUIRE(center[0] == 1);
            REQUIRE(center[1] == 2);
            REQUIRE(center[2] == 3);

            REQUIRE(a.mass() == 1.07);
            REQUIRE(a.charge() == 1);
        }
    }
}

