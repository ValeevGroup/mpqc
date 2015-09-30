#pragma once
#ifndef MPQC_CLUSTERING_FUNCTIONS_H
#define MPQC_CLUSTERING_FUNCTIONS_H

#include "molecule_fwd.h"
#include <vector>

namespace mpqc {
namespace molecule {

namespace clustering {
using output_t = std::vector<Cluster>;
using input_t = std::vector<Clusterable>;

double sum_cluster_distances(const std::vector<Cluster> &clusters);
bool none_zero(std::vector<Cluster> const &clusters);

class kmeans {
  public:
    kmeans(unsigned long seed);
    output_t operator()(input_t const &clusterables, unsigned long nclusters);

  private:
    void initialize_clusters(input_t const &clusterables);

    void attach_clusterables(input_t const &cs);

    std::vector<Vec3D>
    update_clusters(const std::vector<Clusterable> &clusterables);

    bool kmeans_converged(std::vector<Vec3D> const &old_centers);

    output_t cluster(std::vector<Clusterable> const &clusterables);

    std::vector<Cluster>::iterator
    closest_cluster(const std::vector<Cluster>::iterator begin,
                    const std::vector<Cluster>::iterator end,
                    Vec3D const &center);

    unsigned long seed_;
    std::vector<Cluster> clusters_;
};


} // namespace clustering
} // namespace molecule
} // namespace tcc

#endif // MPQC_CLUSTERING_FUNCTIONS_H
