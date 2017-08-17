#include "catch.hpp"

#include <mpqc/chemistry/qc/lcao/wfn/ao_wfn.h>
#include <sstream>

#include "mpqc/chemistry/qc/lcao/factory/lcao_factory.h"

using namespace mpqc;

TEST_CASE("AOWfn", "[aowfn]") {
  const char atoms_xyz_cstr[] =
      "2\n"
      "\n"
      "He 0 0  0\n"
      "Ne 0 0  2\n";

  libint2::initialize();

  using AOWfn = lcao::AOWavefunction<TA::TensorD, TA::SparsePolicy>;
  std::shared_ptr<AOWfn> aowfn;

  SECTION("constructor") {
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

  SECTION("compute ao integrals") {
    REQUIRE_NOTHROW(auto S = aowfn->ao_factory().compute(L"<μ|ν>"));
  }

  libint2::finalize();

}  // Molecule test case
