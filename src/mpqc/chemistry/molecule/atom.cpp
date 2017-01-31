#include "mpqc/chemistry/molecule/atom.h"

#include <iostream>
#include <map>

#include "mpqc/chemistry/units/units.h"
#include <libint2/chemistry/elements.h>

namespace mpqc {

namespace detail {

std::string Z_to_element_name(int64_t Z) {
  using libint2::chemistry::element_info;
  for (const auto &e : element_info) {
    if (e.Z == Z) return e.symbol;
  }
  abort();
}

}  // namespace detail

std::string Atom::xyz_string(bool convert_to_angstroms) const {
  std::string name = detail::Z_to_element_name(atomic_number_);
  name += ' ';

  auto unit_factory = UnitFactory::get_default();
  auto bohr_to_ang = unit_factory->make_unit("angstrom").from_atomic_units();

  Vector3d center = center_;
  if (convert_to_angstroms) {
    center *= bohr_to_ang;
  }

  name += (std::to_string(center[0]) + ' ');
  name += (std::to_string(center[1]) + ' ');
  name += std::to_string(center[2]);

  return name;
}

std::ostream &operator<<(std::ostream &os, Atom const &a) {
  os << "Atom: {Z: " << a.charge() << ", mass: " << a.mass()
     << ", pos: " << a.center().transpose() << "}";
  return os;
}

}  // namespace mpqc
