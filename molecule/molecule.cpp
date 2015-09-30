#include "molecule.h"
#include "atom_based_cluster_concept.h"
#include "atom_based_cluster.h"

#include "common.h"
#include "atom_masses.h"
#include "atom.h"

#include <iostream>

namespace mpqc {
namespace molecule {

namespace molecule_detail {

inline double calculate_mass(const std::vector<Clusterable> &cs) {
    using iter_t = decltype(cs.begin());
    assert(false);  // Fix c.mass() piece
    return 1.0;
   //  return std::accumulate(cs.begin(), cs.end(), 0.0,
   //                         [](double val,
   //                            Clusterable const &c) { return val + c.mass(); });
}

inline int calculate_charge(std::vector<Clusterable> const &cs) {

    assert(false);  // Fix c.mass() piece
    return 1.0;
    // return std::accumulate(cs.begin(), cs.end(), 0.0,
    //                        [](double val,
    //                           Clusterable const &c) { return val + c.charge(); });
}

// Functor for sorting centers based on the distance from a point.
class sort_by_distance_from_point {
  public:
    sort_by_distance_from_point(const Vec3D point) : point_(point) {}

    bool operator()(const Clusterable &a, const Clusterable &b) const {
        Vec3D a_dist = a.center() - point_;
        Vec3D b_dist = b.center() - point_;
        if (!(a_dist.squaredNorm() == b_dist.squaredNorm())) {
            return a_dist.squaredNorm() < b_dist.squaredNorm();
        } else if (a_dist[0] == b_dist[0]) {
            return a_dist[0] < b_dist[0];
        } else if (a_dist[1] == b_dist[1]) {
            return a_dist[1] < b_dist[1];
        } else
            return a_dist[2] < b_dist[2];
    }

  private:
    Vec3D point_;
};

void sort_elements(std::vector<Clusterable> &elems, const Vec3D &point) {
    std::sort(elems.begin(), elems.end(),
                       sort_by_distance_from_point(point));
}
} // namespace moleucle detail


Molecule::Molecule(std::vector<AtomBasedClusterable> c) : elements_(std::move(c)) {
    // mass_ = molecule_detail::calculate_mass(elements_);
    // charge_ = molecule_detail::calculate_charge(elements_);
    // center_ = center_of_mass(elements_, mass_);
    // molecule_detail::sort_elements(elements_, center_);
}

Vec3D const &Molecule::center() const { return center_; }
int64_t Molecule::charge() const { return charge_; }

int64_t Molecule::occupation(unsigned long total_charge) const {
    return charge_ - total_charge;
}

double Molecule::mass() const { return mass_; }

double Molecule::nuclear_repulsion() const {

   //  // Have to get the atoms from each cluster
   //  std::vector<Atom> atoms;
   //  for (auto const &cluster : elements_) {
   //      auto c_atoms = cluster.atoms();
   //      atoms.insert(atoms.end(), c_atoms.begin(), c_atoms.end());
   //  }

   //  double energy = 0.0;
   //  for (auto i = 0ul; i < atoms.size(); ++i) {
   //      for (auto j = i + 1; j < atoms.size(); ++j) {
   //          const auto diff = atoms[i].center() - atoms[j].center();
   //          const auto r = diff.norm();
   //          energy += atoms[i].charge() * atoms[j].charge() / r;
   //      }
   //  }

   //  return energy;
}

std::vector<AtomBasedClusterable>::const_iterator Molecule::begin() const {
    return elements_.begin();
}

std::vector<AtomBasedClusterable>::const_iterator Molecule::end() const {
    return elements_.end();
}

int64_t Molecule::nclusters() const { return elements_.size(); }

int64_t Molecule::core_electrons() const {

 //    // get all the atoms
 //    std::vector<Atom> atoms;
 //    for (auto const &cluster : elements_) {
 //        auto c_atoms = cluster.atoms();
 //        atoms.insert(atoms.end(), c_atoms.begin(), c_atoms.end());
 //    }

 //    // loop over all atoms to get the total core electron number
 //    int i, n = 0;
 //    for (i = 0; i < atoms.size(); i++) {
 //        if (atoms[i].charge() == 0) continue;
 //        int z = atoms[i].charge();
 //        if (z > 2) n += 2;
 //        if (z > 10) n += 8;
 //        if (z > 18) n += 8;
 //        if (z > 30) n += 10;
 //        if (z > 36) n += 8;
 //        if (z > 48) n += 10;
 //        if (z > 54) n += 8;
 //        if (z > 72) {
 //            throw("Molecule::core_electrons: atomic number too large");
 //        }
 //    }
 //    return n;
}// 

Molecule read_xyz(std::string const &file_name) {
 //    std::ifstream xyz_file(file_name);
 //    if (xyz_file.fail()) {
 //        std::ostringstream oss;
 //        oss << "could not open file \"" << file_name << "\"";
 //        throw std::invalid_argument(oss.str().c_str());
 //    }
 //    auto libint_atoms = libint2::read_dotxyz(xyz_file);
 //    xyz_file.close();

 //    std::vector<Clusterable> clusterables;
 //    for (auto const &l_atom : libint_atoms) {
 //        Atom atom({l_atom.x, l_atom.y, l_atom.z},
 //                  masses::masses[l_atom.atomic_number], l_atom.atomic_number);
 //        clusterables.emplace_back(std::move(atom));
 //    }

 //    return Molecule{std::move(clusterables)};
}// 

} // namespace molecule
} // namespace mpqc
