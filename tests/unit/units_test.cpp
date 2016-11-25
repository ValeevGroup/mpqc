//
// Created by Chong Peng on 8/25/16.
//

#include "catch.hpp"
#include "mpqc/chemistry/units/units.h"

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
    auto ang = ufac->make_unit("angstrom");
    auto cm = ufac->make_unit("centimeter");
    REQUIRE(ang.to_atomic_units() == Approx( 1.88973 ));
    //std::cout << "1 angstrom = " << ang.to_atomic_units() << " bohr"<< std::endl;
    REQUIRE(cm.to_atomic_units() == Approx( 1.88973e+08 ));
    //std::cout << "1 cm = " << cm.to_atomic_units() << " bohr"<< std::endl;
    REQUIRE(ang.from_atomic_units() == Approx( 0.529177 ));
    //sstd::cout << "1 bohr = " << ang.from_atomic_units() << " angstrom"<< std::endl;
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
}
