#include "catch.hpp"

#include <mpqc/chemistry/qc/lcao/wfn/ao_wfn.h>
#include <sstream>

#include "mpqc/chemistry/qc/lcao/factory/lcao_factory.h"

using namespace mpqc;

SCENARIO("AOWfn is usable", "[aowfn]") {
  const char atoms_xyz_cstr[] =
      "2\n"
      "\n"
      "He 0 0  0\n"
      "Ne 0 0  2\n";

  libint2::initialize();

  GIVEN("an AOWavefunction") {
    using AOWfn = lcao::AOWavefunction<TA::TensorD, TA::SparsePolicy>;
    std::shared_ptr<AOWfn> aowfn;

    {
      std::stringstream iss((std::string(atoms_xyz_cstr)));
      std::shared_ptr<Molecule> mol;
      REQUIRE_NOTHROW(mol = std::make_shared<Molecule>(iss));

      using AtomicBasis = lcao::gaussian::AtomicBasis;
      std::shared_ptr<AtomicBasis> obs;
      REQUIRE_NOTHROW(obs = std::make_shared<AtomicBasis>(
                          KeyVal()
                              .assign("atoms", mol)
                              .assign("world", &TiledArray::get_default_world())
                              .assign("name", "3-21G")));

      REQUIRE_NOTHROW(aowfn = std::make_shared<AOWfn>(
                          KeyVal()
                              .assign("world", &TiledArray::get_default_world())
                              .assign("atoms", mol)
                              .assign("basis", obs)));
    }

    WHEN("overlap ao integrals are computed") {
      TA::TSpArrayD S;
      REQUIRE_NOTHROW(S = aowfn->ao_factory().compute(L"<μ|ν>"));
      THEN("element {0,0} is 1") {
        auto s_tile_0_0 = S.find({0, 0}).get();
        REQUIRE(s_tile_0_0(0) == Approx(1.0));
      }
    }
  }

  libint2::finalize();

}  // Molecule test case
