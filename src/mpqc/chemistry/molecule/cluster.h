#pragma once
#ifndef MPQC_MOLCULE_CLUSTER_H
#define MPQC_MOLCULE_CLUSTER_H


#include <mpqc/chemistry/molecule/cluster_concept.h>
#include <mpqc/chemistry/molecule/molecule_fwd.h>

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
 * \brief is the unit that holds a collection of clusterables that go together.
 *
 * Cluster will hold a vector of clusterables that all belong together.  
 * To update the center of the cluster compute_com must be called, this is 
 * to avoid computing a new center of mass (COM) every time a clusterable is 
 * added.  
 */
class Cluster {
  private:
    std::vector<Clusterable> elements_;
    Vec3D center_ = {0, 0, 0};

  public:
    Cluster() = default;
    Cluster(const Cluster &c) = default;
    Cluster &operator=(const Cluster &c) = default;

    Cluster(Cluster &&c) = default;
    Cluster &operator=(Cluster &&c) = default;

    Cluster(std::vector<Clusterable> const &elems) : elements_(elems) {}
    Cluster(std::vector<Clusterable> &&elems) : elements_(std::move(elems)) {}

    template<typename ...Cs>
    Cluster(Cs ...cs) : elements_{std::move(cs)...} {}

    template <typename T>
    void add_clusterable(T t) { elements_.emplace_back(std::move(t)); }

    int64_t size() const { return elements_.size(); }

    void clear() { elements_.clear(); }

    /**
     * @brief sets the center equal to a point.
     */
    void set_center(Vec3D point) { center_ = point; }

    /**
     * @brief will update the center based on the current elements.
     */
    void update_center();

    /**
     * @brief sum_distances_from_center calculates the sum of the disances of
     * each clusterable to the center of the cluster.
     * @return reduction over the distances to the cluster center of the
     * clusterables.
     */
    double sum_distances_from_center() const;

    inline Vec3D const& center() const { return center_; }

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
};


// External interface
/*! \brief print the cluster by printing each of its elements
 * 
 */
std::ostream & operator<<(std::ostream &, Cluster const &);

inline Vec3D const& center(Cluster const &c){
    return c.center();
}

/*! @} */

} // namespace molecule
} // namespace mpqc

#endif // MPQC_MOLECULE_CLUSTER_H
