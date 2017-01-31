
#include "catch.hpp"

#include "mpqc/chemistry/molecule/atom.h"

SCENARIO("atoms can be initialized", "[atom]"){
    GIVEN("a default initalized atom"){
        mpqc::Atom a;
        auto center = a.center();
        REQUIRE(center[0] == 0);
        REQUIRE(center[1] == 0);
        REQUIRE(center[2] == 0);

        REQUIRE(a.mass() == 0);
        REQUIRE(a.charge() == 0);

        WHEN("the atom is assigned to a new atom"){
            a = mpqc::Atom({1.0, 2.0, 3.0}, 1.07, 1);
            center = a.center();
            REQUIRE(center[0] == 1);
            REQUIRE(center[1] == 2);
            REQUIRE(center[2] == 3);

            REQUIRE(a.mass() == 1.07);
            REQUIRE(a.charge() == 1);
        }
    }
}

SCENARIO("atoms can be stringified", "[atom]"){
    GIVEN("A hydrogen at 1.0 2.0 3.0"){
        auto a = mpqc::Atom({1.0, 2.0, 3.0}, 1.07, 1);
        std::string s = "H 1.000000 2.000000 3.000000";
        std::string s_bohr = "H 0.529177 1.058354 1.587532";

        // This routine uses std::to_string so formatting is specified by the std.
        REQUIRE(s_bohr == a.xyz_string());
        REQUIRE(s == a.xyz_string(false));
    }
}

