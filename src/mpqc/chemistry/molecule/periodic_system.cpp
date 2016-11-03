#include <mpqc/chemistry/molecule/molecule.h>
#include <mpqc/chemistry/molecule/periodic_system.h>

#include <mpqc/chemistry/molecule/common.h>
#include <mpqc/chemistry/molecule/atom_masses.h>
#include "clustering_functions.h"
#include "mpqc/util/keyval/forcelink.h"

#include <libint2/atom.h>

MPQC_CLASS_EXPORT2("PeriodicSystem", mpqc::molecule::PeriodicSystem);

namespace mpqc {
namespace molecule {

// PeriodicSystem::PeriodicSystem(std::istream &file_stream, Vec3D
// lattice_vector, Vec3I lattice_sum) : Molecule() {
//    init(file_stream, false);
//    dcell_ = lattice_vector;
//    R_max_ = lattice_sum;
//}

PeriodicSystem::PeriodicSystem(const KeyVal &kv) : Molecule(kv) {
    dcell_ = decltype(dcell_)(kv.value<std::vector<double>>("direct_lattice_vector").data());
    const auto angstrom_to_bohr = 1 / 0.52917721092;  // 2010 CODATA value
    dcell_ *= angstrom_to_bohr;
    R_max_ = decltype(R_max_)(kv.value<std::vector<int>>("rmax").data());
    RD_max_ = decltype(RD_max_)(kv.value<std::vector<int>>("rdmax").data());
    RJ_max_ = decltype(RJ_max_)(kv.value<std::vector<int>>("rjmax").data());
    nk_ = decltype(nk_)(kv.value<std::vector<int>>("k_points").data());

    if (R_max_ != RD_max_) {
        throw std::runtime_error("rmax and rdmax do not match! ");
    }
}

double PeriodicSystem::nuclear_repulsion() const {
  auto const &atoms = this->atoms();
  double enuc = 0.0;
  for (auto nx = -RJ_max_(0); nx <= RJ_max_(0); ++nx) {
    for (auto ny = -RJ_max_(1); ny <= RJ_max_(1); ++ny) {
      for (auto nz = -RJ_max_(2); nz <= RJ_max_(2); ++nz) {
        Vec3D shift(nx * dcell_(0),
                    ny * dcell_(1),
                    nz * dcell_(2));

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

void PeriodicSystem::print(std::ostream &os) const {
  os << "Molecule info:" << std::endl;
  os << "\tC.O.M: " << com().transpose() << std::endl;
  os << "\tCharge: " << charge() << std::endl;
  os << "\tMass: " << mass() << std::endl;

  os << "\nElements:\n";
  auto last = end();
  auto second_last = last - 1;
  for (auto it = begin(); it != last; ++it) {
    os << "\t";
    it->print(os) << std::endl;
  }

  os << "\nPeriodic Boundary Conditions:" << std::endl;
  os << "\tDirect lattice vector (in Bohr): [" << dcell_.transpose() << "]"
     << std::endl;
  os << "\tRange of expansion of Bloch Gaussians in AO Gaussians: [" << R_max_.transpose() << "]" << std::endl;
  os << "\tRange of Coulomb operation: [" << RJ_max_.transpose() << "]" << std::endl;
  os << "\tRange of density representation: [" << RD_max_.transpose() << "]" << std::endl;
  os << "\t# of k points in each direction: [" << nk_.transpose() << "]" << std::endl;
}

}  // molecule namespace
}  // mpqc namespace
