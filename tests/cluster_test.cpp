#include "catch.hpp"

#include "../molecule/cluster.h"

#include <sstream>

using namespace tcc;

SCENARIO("Clusters can be printed", "[cluster]") {
    GIVEN("A cluster with 2 hydrogens") {
        auto a = tcc::molecule::Atom({1.0, 2.0, 3.0}, 1.07, 1);
        auto b = tcc::molecule::Atom({-1.0, 2.0, 3.0}, 1.07, 1);

        molecule::Cluster c;
        c.add_clusterable(std::move(a));
        c.add_clusterable(std::move(b));

        std::string s = "2\n\nH 0.529177 1.058354 1.587532\n";
        s+="H -0.529177 1.058354 1.587532\n";


        std::stringstream ss;
        ss << c;

        REQUIRE(ss.str() == s);
    }
}
