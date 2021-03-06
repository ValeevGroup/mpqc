#include "mpqc/chemistry/molecule/common.h"
#include <libint2/atom.h>

#include "mpqc/chemistry/molecule/atom.h"

namespace mpqc {

std::vector<libint2::Atom> to_libint_atom(std::vector<Atom> const &atoms) {
  std::vector<libint2::Atom> l_atoms;
  l_atoms.reserve(atoms.size());

  for (auto const &atom : atoms) {
    const auto &c = atom.center();

    libint2::Atom l_atom;
    l_atom.atomic_number = atom.atomic_number();
    l_atom.x = c[0];
    l_atom.y = c[1];
    l_atom.z = c[2];

    l_atoms.push_back(std::move(l_atom));
  }

  return l_atoms;
}

}  // namespace mpqc
