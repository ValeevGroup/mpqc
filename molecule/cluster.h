#pragma once
#ifndef CLUSTER_H
#define CLUSTER_H

#include <vector>
#include <numeric>
#include "../include/eigen.h"
#include "cluster_concept.h"

namespace tcc {
namespace molecule {
/**
 * @brief The Cluster class contains a vector of clusterables.
 */
class Cluster {
  public:
    using position_t = Clusterable::position_t;

    Cluster() = default;
    Cluster(const Cluster &c) = default;
    Cluster &operator=(const Cluster &c) = default;

    Cluster(Cluster &&c) noexcept;
    Cluster &operator=(Cluster &&c) noexcept;

    /**
     * add_clusterable takes any type which is clusterable and adds it to the
     * cluster.
     */
    template <typename T>
    void add_clusterable(T t) {
        mass_ += t.mass();
        charge_ += t.charge();
        elements_.emplace_back(std::move(t));
    }

    unsigned long nelements() const { return elements_.size(); }

    void clear() {
        charge_ = 0.0;
        mass_ = 0.0;
        elements_.clear();
    }

    /**
     * @brief init_center just sets the center equal to a point.
     * @param point is a vector which will be swapped into the center.
     */
    void init_center(position_t point) { center_.swap(point); }

    /**
     * @brief guess_center will compute a new center of mass for the cluster.
     */
    void guess_center();

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
    inline double charge() const { return charge_; }

    /**
     * @brief begin returns the begin iterator to the vector of clusterables.
     */
    inline std::vector<const Clusterable>::iterator begin() const {
        return elements_.begin();
    }

    /**
     * @brief end returns the end iterator to the vector of clusterables
     */
    inline std::vector<const Clusterable>::iterator end() const {
        return elements_.end();
    }

  private:
    std::vector<Clusterable> elements_;
    position_t center_ = {0, 0, 0};
    double charge_ = 0.0;
    double mass_ = 0.0;
};

} // namespace molecule
} // namespace tcc

#endif // CLUSTER_H
