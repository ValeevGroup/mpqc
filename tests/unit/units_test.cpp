//
// Created by Chong Peng on 8/25/16.
//

#include "catch.hpp"
#include "mpqc/chemistry/units/units.h"
#include "mpqc/util/core/exenv.h"

#include <memory>

TEST_CASE("Units", "[units]") {
  SECTION("UnitsFactory ctors") {
    using mpqc::UnitFactory;
    REQUIRE_NOTHROW(UnitFactory ufac("CODATA2014"));
    REQUIRE_NOTHROW(UnitFactory ufac("2014CODATA"));
    REQUIRE_NOTHROW(UnitFactory ufac("CODATA2010"));
    REQUIRE_NOTHROW(UnitFactory ufac("2010CODATA"));
    REQUIRE_NOTHROW(UnitFactory ufac("MPQC2"));
  }
  SECTION("UnitsFactory singleton") {
    using mpqc::UnitFactory;
    REQUIRE_NOTHROW(UnitFactory::get_default());
    REQUIRE(UnitFactory::get_default()->system() == "2014CODATA");
    REQUIRE_NOTHROW(UnitFactory::set_default("MPQC2"));
    REQUIRE(UnitFactory::get_default()->system() == "MPQC2");
  }
  SECTION("simple units") {
    using mpqc::UnitFactory;
    auto ufac = UnitFactory::get_default();
    auto cm = ufac->make_unit("centimeter");
    REQUIRE(cm.to_atomic_units() == Approx( 1.88973e+08 ));
    //std::cout << "1 cm = " << cm.to_atomic_units() << " bohr"<< std::endl;
    REQUIRE(cm.from_atomic_units() == Approx( 5.29177e-09 ));
    //std::cout << "1 bohr = " << cm.from_atomic_units() << " cm"<< std::endl;
    auto kcal_per_mol = ufac->make_unit("kcal_per_mol");
    REQUIRE(kcal_per_mol.from_atomic_units() == Approx( 627.51 ));
    auto eV = ufac->make_unit("eV");
    REQUIRE(eV.from_atomic_units() == Approx( 27.21138 ));

  }
  SECTION("composite units") {
    using mpqc::UnitFactory;
    auto ufac = UnitFactory::get_default();
    auto m_p_sec = ufac->make_unit("m/s");
    //std::cout << "1 m/s = " << m_p_sec.to_atomic_units() << " a.u."<< std::endl;
    REQUIRE(m_p_sec.to_atomic_units() == Approx( 4.57103e-7 ));
    auto kcal_per_mol = ufac->make_unit("kcal/mol");
    REQUIRE(kcal_per_mol.from_atomic_units() == Approx( 627.51 ));
  }
  SECTION("The Units library synopsis") {
    using namespace mpqc;
    using namespace std;
    // redirect I/O to /dev/null
    auto printnode = FormIO::get_printnode();
    FormIO::set_printnode(-1);
//! [The Units library snippet]
UnitFactory::set_default("2010CODATA");  // will use the 2010 CODATA revision of the physical constants

// simple unit conversion
Unit ang = UnitFactory::get_default()->make_unit("angstrom");
ExEnv::out0() << "1 angstrom = " << ang.to_atomic_units() << " bohr" << endl;
ExEnv::out0() << "1 bohr = " << ang.from_atomic_units() << " angstrom" << endl;

// composite unit conversion
Unit m_p_s2 = UnitFactory::get_default()->make_unit("m / s * s");
ExEnv::out0() << "1 m/s^2 = " << m_p_s2.to_atomic_units() << " a.u. of acceleration" << endl;
ExEnv::out0() << "1 a.u. of acceleration = " << m_p_s2.from_atomic_units() << " m/s^2" << endl;
//! [The Units library snippet]
    // revert I/O
    FormIO::set_printnode(printnode);
    REQUIRE(ang.to_atomic_units() == Approx( 1.88973 ));
    REQUIRE(ang.from_atomic_units() == Approx( 0.529177 ));
    REQUIRE(m_p_s2.to_atomic_units() == Approx( 1.105679e-23 ));
    REQUIRE(m_p_s2.from_atomic_units() == Approx( 9.04422e22 ));
  }
}
