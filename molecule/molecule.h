#pragma once
#ifndef MPQC_MOLECULE_MOLECULE_H
#define MPQC_MOLECULE_MOLECULE_H

#include "molecule_fwd.h"
#include "atom_based_cluster.h"
#include "mpqc/util/keyval/keyval.hpp"

#include <iosfwd>
#include <vector>


namespace mpqc {
namespace molecule {

/*!
 * \defgroup Molecule Molecule
 *
 * \brief The molecule module contains information about how to make and cluster
 *molecules
 *
 *
 *  The molecule module contains all of the classes which are needed for
 *  clustering.  Ultimately, everything in the molecule module should support
 *  the AtomClusterable interface.
 *
 * @{
 */

/*! \brief Molecule is a class which contains a vector of AtomBasedClusterables
 *
 * At its core molecule is a collection of things which can be collapsed
 * to atoms.  Its main job is allow for clustering of its clusterables.
 *
 */
class Molecule : public DescribedClass {
  private:
    std::vector<AtomBasedClusterable> elements_;

    Vec3D com_ = {0, 0, 0};
    double mass_ = 0.0;
    int64_t charge_ = 0;

    void init(std::istream& file, bool sort_input);

  public:

    Molecule() = default;

    Molecule(const KeyVal& kv);

    Molecule(std::istream &file_stream, bool sort_input);

    Molecule(std::vector<AtomBasedClusterable> c, bool sort_input = true);

    void set_mass(double mass) {
        Molecule::mass_ = mass;
    }

    void set_charge(int64_t charge) {
        Molecule::charge_ = charge;
    }

    int64_t charge() const { return charge_; }
    double mass() const { return mass_; }

    std::vector<AtomBasedClusterable>::const_iterator begin() const {
        return elements_.begin();
    }

    std::vector<AtomBasedClusterable>::const_iterator end() const {
        return elements_.end();
    }

    int64_t nclusters() const { return elements_.size(); }

    std::vector<AtomBasedClusterable> const &clusterables() const {
        return elements_;
    }

    int64_t occupation(unsigned long total_charge) const {
        return charge_ - total_charge;
    }

    int64_t core_electrons() const;

    double nuclear_repulsion() const;

    std::vector<Atom> atoms() const;

    Vec3D const &com() const { return com_; }
};

std::ostream &operator<<(std::ostream &, Molecule const &);

/*! @} */

} // namespace molecule
} // namespace mpqc

#endif // MPQC_MOLECULE_MOLECULE_H
