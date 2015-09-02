#include "clustering_functions.h"
#include "common.h"
#include "cluster.h"

#include <cassert>
#include <numeric>
#include <iostream>
#include <random>

namespace tcc {
namespace molecule {
namespace clustering {

double sum_cluster_distances(const std::vector<Cluster> &clusters) {
    return std::accumulate(clusters.begin(), clusters.end(), 0.0,
                           [](double d, const Cluster &c) {
        return d + c.sum_distances_from_center();
    });
}

bool none_zero(std::vector<Cluster> const &clusters) {
    return std::all_of(clusters.begin(), clusters.end(),
                       [](Cluster const &c) { return c.nelements() != 0; });
}

kmeans::kmeans(unsigned long seed) : seed_{seed}, clusters_() {}

output_t kmeans::
operator()(input_t const &clusterables, unsigned long nclusters) {
    assert(clusterables.size() >= nclusters);
    clusters_.resize(nclusters);
    initialize_clusters(clusterables);
    return cluster(clusterables);
}

void kmeans::initialize_clusters(const input_t &clusterables) {

    std::vector<double> weights(clusterables.size(), 1.0);
    std::mt19937 engine(seed_);

    auto end = clusters_.end();
    for (auto it = clusters_.begin(); it != end; ++it) {
        std::discrete_distribution<unsigned int> random_index(weights.begin(),
                                                              weights.end());

        auto idx = random_index(engine);
        auto center_guess = clusterables[idx].center();
        it->set_center(center_guess);

        for(auto i = 0ul; i < clusterables.size(); ++i){
            const auto clusterable_center = clusterables[i].center();

            // Find the closest cluster that has been initialized.
            auto last = it;
            std::advance(last, 1);
            const auto cluster_center
                = closest_cluster(clusters_.begin(), last, clusterable_center)
                      ->center();

            // Calculate weight = dist^2
            weights[i] = diff_squaredNorm(clusterable_center, cluster_center);
        }
    }

    attach_clusterables(clusterables);
}

void kmeans::attach_clusterables(const std::vector<Clusterable> &cs) {
    // Erase the ownership information for each cluster.
    for(auto &cluster : clusters_){
        cluster.clear();
    }

    for (auto const &clusterable : cs) {
        auto closest = closest_cluster(clusters_.begin(), clusters_.end(),
                                       clusterable.center());
        closest->add_clusterable(clusterable);
    }

    for(auto &cluster : clusters_){
        cluster.compute_com();
    }
}

std::vector<Cluster>::iterator
kmeans::closest_cluster(std::vector<Cluster>::iterator const begin,
                        std::vector<Cluster>::iterator const end,
                        position_t const &center) {
    return std::min_element(begin, end,
                            [&](const Cluster &a, const Cluster &b) {
        return diff_squaredNorm(a.center(), center)
               < diff_squaredNorm(b.center(), center);
    });
}


std::vector<position_t>
kmeans::update_clusters(std::vector<Clusterable> const &clusterables) {
    // Vector to hold the previous iterations centers.
    std::vector<position_t> old_centers;
    old_centers.reserve(clusters_.size());

    // Copy the old centers
    for (auto const &cluster : clusters_) {
        old_centers.push_back(cluster.center());
    }

    // Update clusters
    attach_clusterables(clusterables);
    return old_centers;
}

bool kmeans::kmeans_converged(const std::vector<position_t> &old_centers) {
    return std::equal(
        old_centers.begin(), old_centers.end(), clusters_.begin(),
        [](const position_t &old_center, const Cluster &new_cluster) {
            return 1e-7 < (old_center - new_cluster.center()).squaredNorm();
        });
}

output_t kmeans::cluster(const std::vector<Clusterable> &clusterables) {
    unsigned int niter = 100, iter = 0;

    // initialize the old centers and update the clusters
    auto old_centers = update_clusters(clusterables);

    // Compute new clusters and save the old centers for comparison.
    while (!kmeans_converged(old_centers) && (niter > iter++)) {
        old_centers = update_clusters(clusterables);
    }

    return clusters_;
}

} // namespace clustering
} // namespace molecule
} // namespace tcc
