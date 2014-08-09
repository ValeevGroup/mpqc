#ifndef CLUSTERING_FUNCTIONS_H
#define CLUSTERING_FUNCTIONS_H

#include "cluster_concept.h"
#include "cluster.h"
#include <random>

namespace clustering {
using output_t = std::vector<Cluster>;
using input_t = std::vector<Clusterable>;

double sum_cluster_distances(const std::vector<Cluster> &clusters) {
  return std::accumulate(clusters.begin(), clusters.end(), 0.0,
                         [](double d, const Cluster &c) {
    return d + c.sum_distances_from_center();
  });
}

// TODO move this to cpp file
class kmeans {
public:
  kmeans(unsigned long seed) : seed_(std::move(seed)), clusters_() {}
  output_t operator()(input_t clusterables, unsigned long nclusters) {
    clusters_.resize(nclusters);
    initialize_clusters(clusterables);
    return cluster(clusterables);
  }

private:
  void initialize_clusters(const std::vector<Clusterable> &clusterables) {
    std::vector<double> weights(clusterables.size(), 1.0);
    std::default_random_engine engine(seed_);

    for (auto &cluster : clusters_) {
      std::discrete_distribution<unsigned int> random_index(weights.begin(),
                                                            weights.end());

      auto center_guess = clusterables[random_index(engine)].center();
      cluster.init_center(center_guess);

      std::transform(weights.begin(), weights.end(), clusterables.begin(),
                     weights.begin(), [&](double d, const Clusterable &c) {
        double dist = (c.center() - center_guess).norm();
        return (d > 1e-15) ? dist * dist : 0.0;
      });
    }

    attach_clusterables(clusterables);
  }

  void attach_clusterables(const std::vector<Clusterable> &cs) {
    for (auto &cluster : clusters_) {
      cluster.clear();
    }

    for (const auto &c : cs) {
      min_element(clusters_.begin(), clusters_.end(),
                  [&](const Cluster &a, const Cluster &b) {
                    return (a.center() - c.center()).norm() <
                           (b.center() - c.center()).norm();
                  })->add_clusterable(c);
    }
  }

  std::vector<Cluster::position_t>
  update_clusters(const std::vector<Clusterable> &clusterables) {
    std::vector<Clusterable::position_t> old_centers;

    for (auto &cluster : clusters_) {
      old_centers.emplace_back(cluster.center());
      cluster.guess_center();
    }

    attach_clusterables(clusterables);

    return std::move(old_centers);
  }

  bool kmeans_converged(const std::vector<Cluster::position_t> &old_centers) {
    return std::equal(
        old_centers.begin(), old_centers.end(), clusters_.begin(),
        [](const Cluster::position_t &old_center, const Cluster &new_cluster) {
          return 1e-15 < (old_center - new_cluster.center()).norm();
        });
  }

  output_t cluster(const std::vector<Clusterable> &clusterables) {
    unsigned int niter = 100, iter = 0;

    auto old_centers = update_clusters(clusterables);

    while (!kmeans_converged(old_centers) && (niter > iter++)) {
      old_centers = update_clusters(clusterables);
    }

    return clusters_;
  }

  unsigned long seed_;
  std::vector<Cluster> clusters_;
};

} // namespace clustering

#endif // CLUSTERING_FUNCTIONS_H
