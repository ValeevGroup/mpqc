#include "catch.hpp"

#include <sstream>

#include "mpqc/chemistry/molecule/molecule.h"

using namespace mpqc;

TEST_CASE("Molecule", "[molecule]") {

  const char atoms_xyz_cstr[] =
      "6\n"
      "\n"
      "Ne 0 0  0\n"
      "Ne 0 0  2\n"
      "Ne 0 0  4\n"
      "Ne 0 0  6\n"
      "Ne 0 0  8\n"
      "Ne 0 0 10\n";

  SECTION("XYZ construction") {
    Molecule mol;
    std::stringstream iss((std::string(atoms_xyz_cstr)));
    REQUIRE_NOTHROW(mol = Molecule(iss));
  }

  SECTION("atom extraction/update") {
    Molecule mol;
    std::stringstream iss((std::string(atoms_xyz_cstr)));
    REQUIRE_NOTHROW(mol = Molecule(iss));

    std::vector<Atom> atoms = mol.atoms();
    CHECK(atoms.size() == 6);

    atoms[1] = Atom({1., 2., 3.}, 30, 50);
    atoms[5] = Atom({2., 3., 4.}, 40, 60);
    REQUIRE_NOTHROW(mol.update(atoms));

    std::vector<Atom> updated_atoms = mol.atoms();
    CHECK(updated_atoms[1].center() == atoms[1].center());
    CHECK(updated_atoms[5].center() == atoms[5].center());
  }

}  // Molecule test case
