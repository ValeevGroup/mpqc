#include "mpqc/chemistry/molecule/unit_cell.h"

#include <libint2/atom.h>

#include "atomic_data.h"
#include "mpqc/chemistry/molecule/clustering_functions.h"
#include "mpqc/chemistry/molecule/common.h"
#include "mpqc/chemistry/molecule/molecule.h"
#include "mpqc/chemistry/units/units.h"
#include "mpqc/util/keyval/forcelink.h"

namespace mpqc {
UnitCell::UnitCell(const KeyVal &kv) : Molecule(kv) {
  dcell_ =
      decltype(dcell_)(kv.value<std::vector<double>>("lattice_param").data());

  auto unit_factory = UnitFactory::get_default();
  auto angstrom_to_bohr = unit_factory->make_unit("angstrom").to_atomic_units();
  dcell_ *= angstrom_to_bohr;
}

double UnitCell::nuclear_repulsion_energy(Vector3i RJ_max) const {
  auto const &atoms = this->atoms();
  double enuc = 0.0;
  for (auto nx = -RJ_max(0); nx <= RJ_max(0); ++nx) {
    for (auto ny = -RJ_max(1); ny <= RJ_max(1); ++ny) {
      for (auto nz = -RJ_max(2); nz <= RJ_max(2); ++nz) {
        Vector3d shift(nx * dcell_(0), ny * dcell_(1), nz * dcell_(2));

        for (auto i = 0ul; i < atoms.size(); ++i) {
          for (auto j = 0ul; j < atoms.size(); ++j) {
            // exclude self-interactions (intra-cell interactions are still included
            if (!(nx == 0 && ny == 0 && nz == 0 && i == j)) {
              auto r = (atoms[i].center() - atoms[j].center() + shift).norm();
              auto e_cell = atoms[i].charge() * atoms[j].charge() / r;
              enuc += e_cell;
            }
          }
        }
      }
    }
  }

  enuc /= 2;  // takes 2 to tango ... energy per cell is half the interaction energy

  return enuc;
}

std::ostream &operator<<(std::ostream &os, UnitCell const &unitcell) {
    os << "Molecule info:" << std::endl;
    os << "\tC.O.M: " << unitcell.com().transpose() << std::endl;
    os << "\tMass: " << unitcell.mass() << std::endl;

  os << "\nElements:\n";
  auto last = unitcell.end();
  for (auto it = unitcell.begin(); it != last; ++it) {
    os << "\t";
    it->print(os) << std::endl;
  }

  os << "\nUnit cell info:" << std::endl;
  os << "\tLattice parameters (in Bohr): [" << unitcell.dcell().transpose()
     << "]" << std::endl;

  return os;
}

}  // namespace mpqc
