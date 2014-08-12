#include "clustering_functions.h"
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/blocked_range.h>
#include <tbb/spin_mutex.h>

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
  // Erase the ownership information in the cluster.
  tbb::parallel_for_each(clusters_.begin(), clusters_.end(),
                         [](Cluster &c) { c.clear(); });

  // Attach each clusterable to a cluster.
  tbb::spin_mutex myMutex;
  tbb::parallel_for_each(cs.begin(), cs.end(), [&](const Clusterable &c) {
    auto iter = min_element(clusters_.begin(), clusters_.end(),
                            [&](const Cluster &a, const Cluster &b) {
      return (a.center() - c.center()).norm() <
             (b.center() - c.center()).norm();
    });
    tbb::spin_mutex::scoped_lock lock(myMutex);
    iter->add_clusterable(c); // must lock, two could try and add at once.
  });

  // for(auto &cluster : clusters_){
  //  cluster.guess_center();
  //}

  tbb::parallel_for_each(clusters_.begin(), clusters_.end(),
                         [](Cluster &c) { c.guess_center(); });
}

std::vector<Cluster::position_t>
kmeans::update_clusters(const std::vector<Clusterable> &clusterables) {
  // temp to hold the old clusters must allocate for parallel loop.
  std::vector<Clusterable::position_t> old_centers(clusters_.size());

  // store the old centers and guess new ones.
  tbb::parallel_for(tbb::blocked_range<unsigned long>(0, clusters_.size()),
                    [&](const tbb::blocked_range<unsigned long> &r) {
    for (auto i = r.begin(); i != r.end(); ++i) {
      old_centers[i] = clusters_[i].center();
    }
  });

  attach_clusterables(clusterables);

  return old_centers; // Let's go rvo!!;
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
