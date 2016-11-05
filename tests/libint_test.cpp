
#include "catch.hpp"

#include <libint2.hpp>

SCENARIO("libint can be initialized and finalized repeatedly", "[libint]"){
  REQUIRE(libint2::initialized() == false);
  libint2::initialize();
  REQUIRE(libint2::initialized() == true);
  libint2::finalize();
  REQUIRE(libint2::initialized() == false);
  libint2::initialize();
  REQUIRE(libint2::initialized() == true);
  libint2::finalize();
  REQUIRE(libint2::initialized() == false);
}

SCENARIO("libint basis can be constructed", "[libint]"){
  libint2::initialize();
  libint2::Atom h1{1, 0.0, 0.0, 0.0};
  libint2::Atom h2{1, 0.0, 0.0, 1.0};
  std::vector<libint2::Atom> mol = {std::move(h1), std::move(h2)};

  libint2::BasisSet sto3g("STO-3G", mol);

  REQUIRE(sto3g.size() == 2);
}

