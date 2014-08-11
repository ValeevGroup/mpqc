#include "clustering_functions.h"

namespace clustering {

double sum_cluster_distances(const std::vector<Cluster> &clusters) {
  return std::accumulate(clusters.begin(), clusters.end(), 0.0,
                         [](double d, const Cluster &c) {
    return d + c.sum_distances_from_center();
  });
}

kmeans::kmeans(unsigned long seed) : seed_(std::move(seed)), clusters_() {}
output_t kmeans::operator()(input_t clusterables, unsigned long nclusters) {
  clusters_.resize(nclusters);
  initialize_clusters(clusterables);
  return cluster(clusterables);
}

void kmeans::initialize_clusters(const input_t &clusterables) {
  std::vector<double> weights(clusterables.size(), 1.0);
  std::default_random_engine engine(seed_);

  // TODO_PAR tbb this loop
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

  // Go ahead and attach the clusterable to their clusters.
  attach_clusterables(clusterables);
}

void kmeans::attach_clusterables(const std::vector<Clusterable> &cs) {
  // TODO_PAR loop over clusters should also be parallel.
  // Erase the ownership information in the cluster.
  for (auto &cluster : clusters_) {
    cluster.clear();
  }

  // TODO_PAR Make loop over clusterables parallel.
  // Attach each clusterable to a cluster.
  for (const auto &c : cs) {
    min_element(clusters_.begin(), clusters_.end(),
                [&](const Cluster &a, const Cluster &b) {
                  return (a.center() - c.center()).norm() <
                         (b.center() - c.center()).norm();
                })->add_clusterable(c);
  }
}

std::vector<Cluster::position_t>
kmeans::update_clusters(const std::vector<Clusterable> &clusterables) {
  // temp to hold the old clusters
  std::vector<Clusterable::position_t> old_centers;

  //TODO_PAR make this loop parallel.
  // store the old centers and guess new ones.
  for (auto &cluster : clusters_) {
    old_centers.emplace_back(cluster.center());
    cluster.guess_center();
  }

  attach_clusterables(clusterables);

  // Hopeing for rvo
  return old_centers;
}

bool
kmeans::kmeans_converged(const std::vector<Cluster::position_t> &old_centers) {
  return std::equal(
      old_centers.begin(), old_centers.end(), clusters_.begin(),
      [](const Cluster::position_t &old_center, const Cluster &new_cluster) {
        return 1e-15 < (old_center - new_cluster.center()).norm();
      });
}

output_t kmeans::cluster(const std::vector<Clusterable> &clusterables) {
  unsigned int niter = 100, iter = 0;

  // initialize the old centers and update the clusters
  auto old_centers = update_clusters(clusterables);

  while (!kmeans_converged(old_centers) && (niter > iter++)) {
    // save old centers and update the clusters at the same time.
    old_centers = update_clusters(clusterables);
  }

  return clusters_;
}

} // namespace clustering
