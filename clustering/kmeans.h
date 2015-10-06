#pragma once
#ifndef MPQC_CLUSTERING_KMEANS_H
#define MPQC_CLUSTERING_KMEANS_H

#include "../common/typedefs.h"
#include "common.h"

#include <iostream>

#include <vector>
#include <random>

namespace mpqc {
namespace clustering {

class Kmeans {
  private:
    int64_t seed_ = 42;
    int64_t max_iters_ = 100;

  public:
    Kmeans() = default;
    Kmeans(Kmeans const &) = default;
    Kmeans(Kmeans &&) = default;
    Kmeans &operator=(Kmeans const &) = default;
    Kmeans &operator=(Kmeans &&) = default;

    Kmeans(int64_t seed) : seed_(seed) {}
    Kmeans(int64_t seed, int64_t iters) : seed_(seed), max_iters_(iters) {}

    template <typename Cluster, typename Clusterable>
    std::vector<Cluster>
    cluster(std::vector<Clusterable> const &cbls, int64_t nclusters) {
        auto clusters = initial_guess<Cluster>(cbls, nclusters);
        lloyds_algorithm(clusters, cbls);
        return clusters;
    }

  private:
    template <typename Cluster>
    std::vector<Vec3D> centers(std::vector<Cluster> const &clusters) const {
        std::vector<Vec3D> centers;
        centers.reserve(clusters.size());

        for (auto const &cluster : clusters) {
            centers.emplace_back(center(cluster));
        }

        return centers;
    }

    template <typename Cluster, typename Clusterable>
    void attach_clusterables(std::vector<Cluster> &clusters,
                             std::vector<Clusterable> const &cbls) {
        for (auto &cl : clusters) {
            clear(cl);
        }

        for (auto const &cbl : cbls) {
            auto closest = closest_cluster(clusters.begin(), clusters.end(),
                                           center(cbl));
            attach_clusterable(*closest, cbl);
        }

        for (auto &cluster : clusters) {
            update_center(cluster);
        }
    }

    template <typename Cluster, typename Clusterable>
    void lloyds_algorithm(std::vector<Cluster> &clusters,
                          std::vector<Clusterable> const &cbls) {

        // Save old centers
        auto old_centers = centers(clusters);

        // Run an iteration to get new clusters
        attach_clusterables(clusters, cbls);
        auto new_centers = centers(clusters);

        auto iter = 0;
        while (!kmeans_converged(old_centers, new_centers)
               && (max_iters_ > iter++)) {
            old_centers = new_centers;
            attach_clusterables(clusters, cbls);
            new_centers = centers(clusters);
        }
    }

    bool kmeans_converged(std::vector<Vec3D> const &old_centers,
                          std::vector<Vec3D> const &new_centers) {
        return std::equal(old_centers.begin(), old_centers.end(),
                          new_centers.begin(),
                          [](Vec3D const &old_center, Vec3D const &new_center) {
            return 1e-7 < (old_center - new_center).squaredNorm();
        });
    }

    // Init with kmeans++
    template <typename Cluster, typename Clusterable>
    std::vector<Cluster>
    initial_guess(std::vector<Clusterable> const &cbls, int64_t nclusters) {

        std::vector<double> weights(cbls.size(), 1.0);
        std::mt19937 engine(seed_);

        std::vector<Cluster> clusters(nclusters);

        auto end = clusters.end();
        for (auto it = clusters.begin(); it != end; ++it) {

            std::discrete_distribution<unsigned int> random_index(
                  weights.begin(), weights.end());

            // Determine a random center
            auto idx = random_index(engine);
            auto center_guess = center(cbls[idx]);

            // Set the current cluster to that center
            set_center(*it, center_guess);

            auto index = 0;
            for (auto const &clusterable : cbls) {
                const auto current_center = center(clusterable);

                // Find the closest cluster that has been initialized.
                auto first = clusters.begin();
                auto last = it;
                std::advance(last, 1);
                auto test = *closest_cluster(first, last, current_center);
                const auto cloest_cluster_center
                      = center(*closest_cluster(first, last, current_center));

                // Calculate weight = dist^2
                const auto dist = (current_center - cloest_cluster_center);
                weights[index++] = dist.squaredNorm();
            }
        }

        // Updates the clusters with current best guess
        attach_clusterables(clusters, cbls);

        return clusters;
    }
};

} // namespace clustering
} // namespace mpqc

#endif // MPQC_CLUSTERING_KMEANS_H
