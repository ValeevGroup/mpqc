#pragma once
#ifndef TCC_MOLCULE_CLUSTER_H
#define TCC_MOLCULE_CLUSTER_H

#include <vector>
#include <numeric>
#include <iosfwd>

#include "cluster_concept.h"
#include "molecule_fwd.h"

namespace tcc {
namespace molecule {
/*! 
 * \ingroup Molecule
 *
 * @{
 */

/*!
 * \brief is the unit that holds a collection of clusterables that go together.
 *
 * Cluster will hold a vector of clusterables that all belong together.  
 * To update the center of the cluster compute_com must be called, this is 
 * to avoid computing a new center of mass (COM) every time a clusterable is 
 * added.  
 */
class Cluster {
  public:
    Cluster() = default;
    Cluster(const Cluster &c) = default;
    Cluster &operator=(const Cluster &c) = default;

    Cluster(Cluster &&c) = default;
    Cluster &operator=(Cluster &&c) = default;

    template <typename T>
    void add_clusterable(T t) {
        mass_ += t.mass();
        charge_ += t.charge();
        elements_.emplace_back(std::move(t));
    }

    unsigned long nelements() const { return elements_.size(); }

    void clear() {
        charge_ = 0;
        mass_ = 0.0;
        elements_.clear();
    }

    /**
     * @brief init_center just sets the center equal to a point.
     * @param point is a vector which will be swapped into the center.
     */
    void set_center(position_t point) { center_.swap(point); }

    /**
     * @brief guess_center will compute a new center of mass for the cluster.
     */
    void compute_com();

    /**
     * @brief sum_distances_from_center calculates the sum of the disances of
     * each
     * clusterable to the center of the cluster.
     * @return reduction over the distances to the cluster center of the
     * clusterables.
     */
    double sum_distances_from_center() const;

    inline position_t center() const { return center_; }
    inline double mass() const { return mass_; }
    inline int charge() const { return charge_; }

    /**
     * @brief begin returns the begin iterator to the vector of clusterables.
     */
    inline std::vector<Clusterable>::const_iterator begin() const {
        return elements_.begin();
    }

    /**
     * @brief end returns the end iterator to the vector of clusterables
     */
    inline std::vector<Clusterable>::const_iterator end() const {
        return elements_.end();
    }

  private:
    std::vector<Clusterable> elements_;
    position_t center_ = {0, 0, 0};
    int charge_ = 0;
    double mass_ = 0.0;
};

/*! \brief print the cluster by collapsing it to atoms
 * 
 * This function prints the cluster by calling collapse to atoms and then 
 * printing the atoms. The clusters are printed in xyz format.
 */
std::ostream & operator<<(std::ostream &, Cluster const &);

/*! @} */

} // namespace molecule
} // namespace tcc

#endif // TCC_MOLECULE_CLUSTER_H
