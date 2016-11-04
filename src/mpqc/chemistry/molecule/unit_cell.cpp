#include "mpqc/chemistry/molecule/unit_cell.h"

#include <libint2/atom.h>

#include "mpqc/chemistry/molecule/molecule.h"
#include "mpqc/chemistry/molecule/atom_masses.h"
#include "mpqc/chemistry/molecule/common.h"
#include "clustering_functions.h"
#include "mpqc/util/keyval/forcelink.h"


MPQC_CLASS_EXPORT2("UnitCell", mpqc::UnitCell);

namespace mpqc {
UnitCell::UnitCell(const KeyVal &kv) : Molecule(kv) {
  dcell_ = decltype(dcell_)(
      kv.value<std::vector<double>>("lattice_param").data());
  const auto angstrom_to_bohr = 1 / 0.52917721092;  // 2010 CODATA value
  dcell_ *= angstrom_to_bohr;
}

double UnitCell::nuclear_repulsion(Vector3i RJ_max) const {
  auto const &atoms = this->atoms();
  double enuc = 0.0;
  for (auto nx = -RJ_max(0); nx <= RJ_max(0); ++nx) {
    for (auto ny = -RJ_max(1); ny <= RJ_max(1); ++ny) {
      for (auto nz = -RJ_max(2); nz <= RJ_max(2); ++nz) {
        Vector3d shift(nx * dcell_(0), ny * dcell_(1), nz * dcell_(2));

        for (auto i = 0ul; i < atoms.size(); ++i) {
          for (auto j = 0ul; j < atoms.size(); ++j) {
            if (nx == 0 && ny == 0 && nz == 0 && i == j)
              enuc += 0.0;
            else {
              auto r = (atoms[i].center() - atoms[j].center() + shift).norm();
              auto e_cell = 0.5 * atoms[i].charge() * atoms[j].charge() / r;
              enuc += e_cell;
            }
          }
        }
      }
    }
  }

  return enuc;
}

void UnitCell::print(std::ostream &os) const {
  os << "Molecule info:" << std::endl;
  os << "\tC.O.M: " << com().transpose() << std::endl;
  os << "\tCharge: " << charge() << std::endl;
  os << "\tMass: " << mass() << std::endl;

  os << "\nElements:\n";
  auto last = end();
  for (auto it = begin(); it != last; ++it) {
    os << "\t";
    it->print(os) << std::endl;
  }

  os << "\nUnit cell info:" << std::endl;
  os << "\tLattice parameters (in Bohr): [" << dcell_.transpose() << "]"
     << std::endl;
}
}  // mpqc namespace
