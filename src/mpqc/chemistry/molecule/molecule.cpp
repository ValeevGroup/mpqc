#include "mpqc/chemistry/molecule/molecule.h"
//#include <madness/world/world.h>

#include "clustering_functions.h"
#include "mpqc/chemistry/molecule/atom_masses.h"
#include "mpqc/chemistry/molecule/common.h"
#include "mpqc/util/external/madworld/parallel_file.h"

#include <libint2/atom.h>

MPQC_CLASS_EXPORT2("Molecule", mpqc::Molecule);

namespace mpqc {

namespace {

using ABCbl = AtomBasedClusterable;

// Functor for sorting centers based on the distance from a point.
class sort_by_distance_from_point {
 public:
  sort_by_distance_from_point(Vector3d const &point) : point_(point) {}

  bool operator()(ABCbl const &a, ABCbl const &b) const {
    Vector3d a_diff = center_of_mass(a) - point_;
    Vector3d b_diff = center_of_mass(b) - point_;

    const auto a_dist2 = a_diff.squaredNorm();
    const auto b_dist2 = b_diff.squaredNorm();

    if (a_dist2 != b_dist2) {  // Not same distance
      return a_dist2 < b_dist2;
    }

    // a and b were the same distance from point
    // so we need to check each dim separately
    // check in x then y then z order
    if (a_diff[0] != b_diff[0]) {
      return a_diff[0] < b_diff[0];
    } else if (a_diff[1] != b_diff[1]) {
      return a_diff[1] < b_diff[1];
    } else {
      return a_diff[2] < b_diff[2];
    }
  }

 private:
  Vector3d point_;
};

void sort_elements(std::vector<ABCbl> &elems, const Vector3d &point) {
  std::sort(elems.begin(), elems.end(), sort_by_distance_from_point(point));
}

}  // namespace anonymous

Molecule::Molecule(std::vector<ABCbl> c, bool sort_input)
    : elements_(std::move(c)),
      com_(center_of_mass(AtomBasedCluster(elements_))),
      mass_(molecule::sum_mass(elements_)),
      charge_(0),
      total_charge_(molecule::sum_charge(elements_)) {
  if (sort_input) {
    sort_elements(elements_, com_);
  }
}

Molecule::Molecule(const KeyVal &kv) {
  //  std::cout << "Construct Molecule" << std::endl;
  auto file_name = kv.value<std::string>("file_name", "");

  // find world one level higher
  madness::World *world = kv.value<madness::World *>("$:world");

  std::stringstream file;
  utility::parallel_read_file(*world, file_name, file);

  bool sort_input = kv.value<bool>("sort_input", true);
  bool sort_origin = kv.value<bool>("sort_origin", false);

  if (sort_origin) {
    init(file, {0.0, 0.0, 0.0});
  } else {
    init(file, sort_input);
  }

  int n_cluster = kv.value<int>("n_cluster", 0);
  bool attach_hydrogen = kv.value<bool>("attach_hydrogen", true);
  // cluster molecule
  if (n_cluster != 0) {
    Molecule clustered_mol;
    if (attach_hydrogen) {
      clustered_mol = attach_hydrogens_and_kmeans(clusterables(), n_cluster);
    } else {
      clustered_mol = kmeans(clusterables(), n_cluster);
    }

    elements_ = std::move(clustered_mol.elements_);
    com_ = std::move(clustered_mol.com_);
    mass_ = std::move(clustered_mol.mass_);
    total_charge_ = std::move(clustered_mol.total_charge_);
  }

  // attention, has to get charge at the end
  charge_ = kv.value<int>("charge", 0);
  if (charge_ > total_charge_) {
    throw std::invalid_argument("Charge > Total Charge of Molecule! \n");
  }
}

Molecule::Molecule(std::istream &file_stream, bool sort_input) {
  init(file_stream, sort_input);
}

Molecule::Molecule(std::istream &file_stream, Vector3d const &point) {
  init(file_stream, point);
}

Molecule::~Molecule() = default;

void Molecule::init(std::istream &file, bool sort_input) {
  auto libint_atoms = libint2::read_dotxyz(file);

  using ABCbl = AtomBasedClusterable;
  std::vector<ABCbl> atoms;
  for (auto const &l_atom : libint_atoms) {
    Atom atom({l_atom.x, l_atom.y, l_atom.z},
              molecule::masses::masses[l_atom.atomic_number],
              l_atom.atomic_number);
    atoms.emplace_back(std::move(atom));
  }

  elements_ = std::move(atoms);
  com_ = molecule::center_of_mass(elements_);
  mass_ = molecule::sum_mass(elements_);
  total_charge_ = molecule::sum_charge(elements_);

  if (sort_input) {
    sort_elements(elements_, com_);
  }
}

void Molecule::init(std::istream &file, Vector3d const &point) {
  auto libint_atoms = libint2::read_dotxyz(file);

  using ABCbl = AtomBasedClusterable;
  std::vector<ABCbl> atoms;
  for (auto const &l_atom : libint_atoms) {
    Atom atom({l_atom.x, l_atom.y, l_atom.z},
              molecule::masses::masses[l_atom.atomic_number],
              l_atom.atomic_number);
    atoms.emplace_back(std::move(atom));
  }

  elements_ = std::move(atoms);
  com_ = molecule::center_of_mass(elements_);
  mass_ = molecule::sum_mass(elements_);
  total_charge_ = molecule::sum_charge(elements_);

  sort_elements(elements_, point);
}

void Molecule::sort_from_point(Vector3d const &point) {
  sort_elements(elements_, point);
}

std::vector<Atom> Molecule::atoms() const {
  return collapse_to_atoms(elements_);
}

double Molecule::nuclear_repulsion() const {
  auto const &atoms = this->atoms();

  double energy = 0.0;
  for (auto i = 0ul; i < atoms.size(); ++i) {
    for (auto j = i + 1; j < atoms.size(); ++j) {
      const auto r = (atoms[i].center() - atoms[j].center()).norm();
      energy += (atoms[i].charge() * atoms[j].charge()) / r;
    }
  }

  return energy;
}

int64_t Molecule::core_electrons() const {
  int64_t n = 0;
  for (auto const &a : this->atoms()) {
    int z = a.charge();
    assert(z != 0);

    if (z > 2) n += 2;
    if (z > 10) n += 8;
    if (z > 18) n += 8;
    if (z > 30) n += 10;
    if (z > 36) n += 8;
    if (z > 48) n += 10;
    if (z > 54) n += 8;
    if (z > 72) {
      throw("Molecule::core_electrons: atomic number too large");
    }
  }
  return n;
}

void Molecule::print(std::ostream & os) const {
    os << "Molecule C.O.M: " << com().transpose() << ", ";
    os << "charge: " << charge() << ", ";
    os << "mass: " << mass() << ", with Elements: {";

    auto last = end();
    auto second_last = last - 1;
    for (auto it = begin(); it != last; ++it) {
      if (it != second_last) {
        it->print(os) << ", ";
      } else {
        it->print(os) << "}";
      }
    }
}

std::ostream &operator<<(std::ostream &os, Molecule const &mol) {
  mol.print(os);
  return os;
}

}  // namespace mpqc
