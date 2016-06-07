#include <mpqc/chemistry/molecule/molecule.h>
//#include <madness/world/world.h>

#include <mpqc/chemistry/molecule/common.h>
#include <mpqc/chemistry/molecule/atom_masses.h>
#include "../../../../utility/parallel_file.h"

#include <libint2/atom.h>


MPQC_CLASS_EXPORT_KEY2(mpqc::molecule::Molecule, "Molecule");

namespace mpqc {
namespace molecule {

namespace {

using ABCbl = AtomBasedClusterable;

// Functor for sorting centers based on the distance from a point.
class sort_by_distance_from_point {
  public:
    sort_by_distance_from_point(Vec3D const &point) : point_(point) {}

    bool operator()(ABCbl const &a, ABCbl const &b) const {

        Vec3D a_diff = center_of_mass(a) - point_;
        Vec3D b_diff = center_of_mass(b) - point_;

        const auto a_dist2 = a_diff.squaredNorm();
        const auto b_dist2 = b_diff.squaredNorm();


        if (a_dist2 != b_dist2) { // Not same distance
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
    Vec3D point_;
};

void sort_elements(std::vector<ABCbl> &elems, const Vec3D &point) {
    std::sort(elems.begin(), elems.end(), sort_by_distance_from_point(point));
}

} // namespace anonymous


Molecule::Molecule(std::vector<ABCbl> c, bool sort_input)
        : elements_(std::move(c)),
          com_(center_of_mass(elements_)),
          mass_(sum_mass(elements_)),
          charge_(sum_charge(elements_)) {
    if (sort_input) {
        sort_elements(elements_, com_);
    }
}

Molecule::Molecule(const KeyVal& kv) {
//  std::cout << "Construct Molecule" << std::endl;
  auto file_name = kv.value<std::string>("file_name", "");

  // find world one level higher
  madness::World* world = kv.value<madness::World*>("..:world");

  std::stringstream file = utility::parallel_read_file(*world,file_name);

  auto sort_input = kv.value<bool>("sort_input", true);


  init(file, sort_input);
}

Molecule::Molecule(std::istream &file_stream, bool sort_input) {
  init(file_stream, sort_input);
}

void
Molecule::init(std::istream& file, bool sort_input) {
  auto libint_atoms = libint2::read_dotxyz(file);

  using ABCbl = AtomBasedClusterable;
  std::vector<ABCbl> atoms;
  for (auto const &l_atom : libint_atoms) {
      Atom atom({l_atom.x, l_atom.y, l_atom.z},
                masses::masses[l_atom.atomic_number], l_atom.atomic_number);
      atoms.emplace_back(std::move(atom));
  }

  elements_ = std::move(atoms);
  com_ = center_of_mass(elements_);
  mass_ = sum_mass(elements_);
  charge_ = sum_charge(elements_);

  if (sort_input) {
    sort_elements(elements_, com_);
  }
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

std::ostream &operator<<(std::ostream &os, Molecule const &mol) {
    os << "Molecule C.O.M: " << mol.com().transpose() << ", ";
    os << "charge: " << mol.charge() << ", ";
    os << "mass: " << mol.mass() << ", with Elements: {";

    auto last = mol.end();
    auto second_last = last - 1;
    for (auto it = mol.begin(); it != last; ++it) {
        if (it != second_last) {
            it->print(os) << ", ";
        } else {
            it->print(os) << "}";
        }
    }

    return os;
}

} // namespace molecule
} // namespace mpqc
