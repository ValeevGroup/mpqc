#include "molecule.h"
#include "cluster_concept.h"
#include "clustering_functions.h"
#include "cluster.h"
#include "attach_hydrogens.h"
#include "../include/tbb.h"
#include "../include/libint.h"
#include "common.h"
#include "atom_masses.h"

#include <functional>
#include <limits>
#include <iostream>
#include <iostream>

namespace tcc {
namespace molecule {

namespace molecule_detail {

inline double calculate_mass(const std::vector<Clusterable> &cs) {
    using iter_t = decltype(cs.begin());
    return tbb::parallel_reduce(
        tbb::blocked_range<iter_t>(cs.begin(), cs.end()), 0.0,
        [](const tbb::blocked_range<iter_t> &r, double d) -> double {
            return std::accumulate(
                r.begin(), r.end(), d,
                [](double x, const Clusterable &c) { return x + c.mass(); });
        },
        std::plus<double>());
}

inline int calculate_charge(std::vector<Clusterable> const &cs) {

    using iter_t = decltype(cs.begin());
    return tbb::parallel_reduce(
        tbb::blocked_range<iter_t>(cs.begin(), cs.end()), 0,
        [](const tbb::blocked_range<iter_t> &r, int d) -> int {
            return std::accumulate(
                r.begin(), r.end(), d,
                [](int x, const Clusterable &c) { return x + c.charge(); });
        },
        std::plus<int>());
}

// Functor for sorting centers based on the distance from a point.
class sort_by_distance_from_point {
  public:
    sort_by_distance_from_point(const position_t point) : point_(point) {}

    bool operator()(const Clusterable &a, const Clusterable &b) const {
        position_t a_dist = a.center() - point_;
        position_t b_dist = b.center() - point_;
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
    position_t point_;
};

void sort_elements(std::vector<Clusterable> &elems, const position_t &point) {
    tbb::parallel_sort(elems.begin(), elems.end(),
                       sort_by_distance_from_point(point));
}
} // namespace moleucle detail


Molecule::Molecule(std::vector<Clusterable> c) : elements_(std::move(c)) {
    mass_ = molecule_detail::calculate_mass(elements_);
    charge_ = molecule_detail::calculate_charge(elements_);
    center_ = center_of_mass(elements_, mass_);
    molecule_detail::sort_elements(elements_, center_);
}

position_t Molecule::center() const { return center_; }
int Molecule::charge() const { return charge_; }
unsigned long Molecule::occupation(unsigned long total_charge) const {
    return charge_ - total_charge;
}
double Molecule::mass() const { return mass_; }

double Molecule::nuclear_repulsion() const {

    // Have to get the atoms from each cluster
    std::vector<Atom> atoms;
    for(auto const &cluster : elements_){
        auto c_atoms = cluster.atoms();
        atoms.insert(atoms.end(), c_atoms.begin(), c_atoms.end());
    }

    double energy = 0.0;
    for(auto i = 0ul; i < atoms.size(); ++i){
        for(auto j = i + 1; j < atoms.size(); ++j){
            const auto diff = atoms[i].center() - atoms[j].center();
            const auto r = diff.norm();
            energy += atoms[i].charge() * atoms[j].charge() / r;
        }
    }

    return energy;
}

std::vector<Clusterable>::const_iterator Molecule::begin() const {
    return elements_.begin();
}

std::vector<Clusterable>::const_iterator Molecule::end() const {
    return elements_.end();
}

unsigned long Molecule::nelements() const { return elements_.size(); }

std::vector<Cluster>
Molecule::cluster_molecule(cluster_fn_t fn, unsigned long nclusters) const {
    return fn(elements_, nclusters);
}

std::vector<Cluster> Molecule::attach_hydrogens() const {
    return clustering::attach_hydrogens()(elements_);
}

std::vector<Cluster>
Molecule::attach_H_and_kmeans(unsigned long nclusters,
                              unsigned long init_seed) const {
    // if we asked for more clusters than possible return the maximum.
    auto h_attached_clusters = attach_hydrogens();
    if (nclusters >= h_attached_clusters.size()) {
        return h_attached_clusters;
    }

    std::vector<Clusterable> new_clusterables;
    new_clusterables.reserve(h_attached_clusters.size());
    for (auto &&c : h_attached_clusters) {
        new_clusterables.push_back(std::move(c));
    }

    // Compute k-means with a new seed and if that seed give a better answer
    // store the seed
    std::vector<Cluster> clusters;
    int best_seed = init_seed;
    double smallest_dist = std::numeric_limits<double>::max();
    const auto nguesses = 10;
    for (auto i = 0; i < nguesses; ++i) {
        clustering::kmeans func(init_seed);
        clusters = func(new_clusterables, nclusters);

        if (clustering::none_zero(clusters)) {
            double new_dist = clustering::sum_cluster_distances(clusters);
            if (new_dist < smallest_dist) {
                smallest_dist = new_dist;
                best_seed = init_seed;
            }
        }
        init_seed += 10 * (i + 1); // Vary the seed by a resonably large amount.
    }

    // Use the best seed to compute the clusters.
    clustering::kmeans func(best_seed);
    return func(new_clusterables, nclusters);
}

Molecule read_xyz(std::string const &file_name) {
    std::ifstream xyz_file(file_name);
    auto libint_atoms = libint2::read_dotxyz(xyz_file);
    xyz_file.close();

    std::vector<Clusterable> clusterables;
    for (auto const &l_atom : libint_atoms) {
        Atom atom({l_atom.x, l_atom.y, l_atom.z},
                  masses::masses[l_atom.atomic_number], l_atom.atomic_number);
        clusterables.emplace_back(std::move(atom));
    }

    return Molecule{std::move(clusterables)};
}

} // namespace molecule
} // namespace tcc