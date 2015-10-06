#pragma once
#ifndef MPQC_MOLCULE_ATOMBASEDCLUSTER_H
#define MPQC_MOLCULE_ATOMBASEDCLUSTER_H

#include "atom_based_cluster_concept.h"
#include "molecule_fwd.h"

#include <vector>
#include <iosfwd>

namespace mpqc {
namespace molecule {
/*!
 * \ingroup Molecule
 *
 * @{
 */

/*!
 * \brief is the unit that holds a collection of atom based clusterables that go
 *together.
 *
 * AtomBasedCluster will hold a vector of clusterables that all belong together.
 * To update the center of the cluster compute_com must be called, this is
 * to avoid computing a new center of mass (COM) every time a clusterable is
 * added.
 */
class AtomBasedCluster {
  private:
    std::vector<AtomBasedClusterable> elements_;
    Vec3D com_ = {0, 0, 0};
    double mass_ = 0.0;
    int64_t charge_ = 0.0;

  public:
    AtomBasedCluster() = default;
    AtomBasedCluster(const AtomBasedCluster &c) = default;
    AtomBasedCluster &operator=(const AtomBasedCluster &c) = default;

    AtomBasedCluster(AtomBasedCluster &&c) = default;
    AtomBasedCluster &operator=(AtomBasedCluster &&c) = default;

    AtomBasedCluster(std::vector<AtomBasedClusterable> const &elems)
            : elements_(elems) {}
    AtomBasedCluster(std::vector<AtomBasedClusterable> &&elems)
            : elements_(std::move(elems)) {}

    // When constructed from list update immediately
    template <typename... Cs>
    AtomBasedCluster(Cs... cs)
            : elements_{std::move(cs)...} {
        update_cluster();
    }

    template <typename T>
    void add_clusterable(T t) {
        elements_.emplace_back(std::move(t));
    }

    int64_t size() const { return elements_.size(); }

    int64_t charge() const { return charge_; }
    double mass() const { return mass_; }

    std::vector<Atom> atoms() const;

    void clear() { elements_.clear(); }

    /**
     * @brief sets the center equal to a point.
     */
    void set_com(Vec3D point) { com_ = point; }

    /**
     * @brief will update the center based on the current elements.
     *
     * This is done as a separate step because it would be inefficent to update
     *after each addition of a AtomBasedClusterable
     */
    void update_cluster();

    inline Vec3D const &com() const { return com_; }

    /**
     * @brief begin returns the begin iterator to the vector of clusterables.
     */
    inline std::vector<AtomBasedClusterable>::const_iterator begin() const {
        return elements_.begin();
    }

    /**
     * @brief end returns the end iterator to the vector of clusterables
     */
    inline std::vector<AtomBasedClusterable>::const_iterator end() const {
        return elements_.end();
    }
};

// External interface
/*! \brief print the cluster by printing each of its elements
 *
 */
std::ostream &operator<<(std::ostream &, AtomBasedCluster const &);

/*! \brief Returns the Center of mass of the cluster.
 *
 * This function exists to allow interfacing with generic clustering code.
 * Overloading center to return the center of mass is in some sense specializing
 * center for atoms.
 */
inline Vec3D const &center(AtomBasedCluster const &c) { return c.com(); }

inline double mass(AtomBasedCluster const &c) { return c.mass(); }

inline int64_t charge(AtomBasedCluster const &c) { return c.charge(); }

inline Vec3D const &center_of_mass(AtomBasedCluster const &c) {
    return c.com();
}

std::vector<Atom> collapse_to_atoms(AtomBasedCluster const&);

inline void set_center(AtomBasedCluster &c, Vec3D const &point){
    c.set_com(point);
}

inline void clear(AtomBasedCluster &c){
    c.clear();
}

template <typename T>
inline void attach_clusterable(AtomBasedCluster &c, T t){
    c.add_clusterable(std::move(t));
}

inline void update_center(AtomBasedCluster &c){
    c.update_cluster();
}


/*! @} */

} // namespace molecule
} // namespace mpqc

#endif // MPQC_MOLECULE_ATOMBASEDCLUSTER_H
