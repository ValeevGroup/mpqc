
#ifndef MPQC_CLUSTERING_KMEANS_H
#define MPQC_CLUSTERING_KMEANS_H

#include "mpqc/math/clustering/common.h"

#include <iostream>

#include <random>
#include <vector>

namespace mpqc {
namespace clustering {
/*!
 * \defgroup Clustering Clustering
 *
 * \brief The clustering module contains classes and functions for generic
 *     clustering
 *
 *  Clustering is the process of taking a collection of things and determining
 *  how to group them in a logical way.
 *
 * @{
 */

/*! \brief Kmeans is a class that performs k-means clustering
 *
 * Kmeans uses k-means++ for initialization of the clusters and then used
 * Lloyd's algorithm to refine the clusters. Kmeans is considered converged
 * when the positions of the centers stop changing or
 * the maximum number of iterations was reached, during Lloyd's algorithm.
 *
 */
class Kmeans {
 private:
  int64_t seed_ = 42;        //!< Seed for RNG used in k-means++.
  int64_t max_iters_ = 100;  //!< Max iterations in Lloyd's algorithm.

 public:
  Kmeans() = default;
  Kmeans(Kmeans const &) = default;
  Kmeans(Kmeans &&) = default;
  Kmeans &operator=(Kmeans const &) = default;
  Kmeans &operator=(Kmeans &&) = default;

  /*! Constructor to Change the initial seed for Kmeans++
   *
   * \param seed int64_t value to set the seed used in k-means++
   */
  Kmeans(int64_t seed) : seed_(seed) {}

  /*!
   * Constructor to change the initial seed for Kmeans++ and the maximum
   * number of iterations used in clustering
   *
   * \param seed int64_t value to set the seed used in k-means++
   * \param iters int64_t value to set the maximum number of iterations used
   *in
   *      Lloyd's algorithm.
   */
  Kmeans(int64_t seed, int64_t iters) : seed_(seed), max_iters_(iters) {}

  /*! \brief Computes Clusters from Clusterables
   *
   *  \note may throw if the number of clusters
   *      does not make sense, i.e. is too small (<= 0) or
   *      too large (> # input Clusterables)
   *
   *  \param cbls A vector of Clusterables to be clustered. Requirements are
   *      listed below
   *
   *  \param nclusters int64_t which specifies how many clusters are to be
   *      computed.
   *
   * Requirements on Clusterable:
   *  - non-intrusive center: returns a Vector3d
   *
   * Requirements on Cluster:
   *  - non-intrusive center: returns a Vector3d.
   *
   *  - non-intrusive set_center: takes Vector3d allowing for the cluster's
   *      center to be set manually, this is needed for initialization
   *      of the clusters with a guess.
   *
   *  - non-intrusive remove_clusterables: removes the clusterables in the
   *      cluster, but leaves the center and anything else, that is needed for
   *      determining distance to clusterables intact.
   *
   *  - non-intrusive attach_clusterable: which adds a new clusterable to the
   *      cluster.  attach_clusterable should not update the center of the
   *      cluster since that must be done in a separate step in Kmeans
   *
   *  - non-intrusive update_center: Should use the cluster's current
   *      clusterables to update the center of the cluster, but should not
   *      remove any clusterables.
   *
   *  - begin and end functions which allow iteration over the elements of the
   *      cluster.
   */
  template <typename Cluster, typename Clusterable>
  std::vector<Cluster> cluster(std::vector<Clusterable> const &cbls,
                               int64_t nclusters) {
    if (nclusters <= 0 || static_cast<std::size_t>(nclusters) > cbls.size()) {
      throw std::invalid_argument(
          "Number of clusters in Kmeans must be "
          "between 1 and the number of elements "
          "to be clustered.");
    }

    auto clusters = initial_guess<Cluster>(cbls, nclusters);
    lloyds_algorithm(clusters, cbls);
    return clusters;
  }

 private:
  // Save old centers so we can see if Kmeans has converged.
  template <typename Cluster>
  std::vector<Vector3d> centers(std::vector<Cluster> const &clusters) const {
    std::vector<Vector3d> centers;
    centers.reserve(clusters.size());

    for (auto const &cluster : clusters) {
      centers.emplace_back(center(cluster));
    }

    return centers;
  }

  template <typename Cluster, typename Clusterable>
  void attach_clusterables(std::vector<Cluster> &clusters,
                           std::vector<Clusterable> const &cbls) {
    // Clear the clusters of their clusterables for the next iteration of
    // Lloyd's algorithm
    for (auto &c : clusters) {
      remove_clusterables(c);
    }

    // Attach each clusterable to the closest cluster.
    for (auto const &cbl : cbls) {
      auto closest =
          closest_cluster(clusters.begin(), clusters.end(), center(cbl));
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
    while (!kmeans_converged(old_centers, new_centers) &&
           (max_iters_ > iter++)) {
      old_centers = new_centers;
      attach_clusterables(clusters, cbls);
      new_centers = centers(clusters);
    }
  }

  bool kmeans_converged(std::vector<Vector3d> const &old_centers,
                        std::vector<Vector3d> const &new_centers) {
    auto compare_center_position = [](Vector3d const &oldc,
                                      Vector3d const &newc) {
      return 1e-7 < (oldc - newc).squaredNorm();
    };

    return std::equal(old_centers.begin(), old_centers.end(),
                      new_centers.begin(), compare_center_position);
  }

  // Init with kmeans++ look up paper
  //  k-means++: The Advantages of Careful Seeding,
  //  by: David Arthur and Sergei Vassilvitskii
  template <typename Cluster, typename Clusterable>
  std::vector<Cluster> initial_guess(std::vector<Clusterable> const &cbls,
                                     int64_t nclusters) {
    // initially each center has the same weight
    std::vector<double> weights(cbls.size(), 1.0);
    std::mt19937 engine(seed_);

    std::vector<Cluster> clusters(nclusters);

    auto end = clusters.end();
    for (auto it = clusters.begin(); it != end; ++it) {
      std::discrete_distribution<unsigned int> random_index(weights.begin(),
                                                            weights.end());

      // Determine a random center
      auto idx = random_index(engine);
      auto center_guess = center(cbls[idx]);

      // Set the current cluster to that center
      set_center(*it, center_guess);

      // Update the weights based on the chosen cluster guesses
      auto last_init_cluster = it;
      std::advance(last_init_cluster, 1);
      update_weights(clusters.begin(), last_init_cluster, cbls, weights);
    }

    // Updates the clusters with current best guess
    attach_clusterables(clusters, cbls);

    return clusters;
  }

  // Function to update the weights in the k-means++ initialization.
  template <typename Iter, typename Clusterable>
  void update_weights(Iter first, Iter last,
                      std::vector<Clusterable> const &cbls,
                      std::vector<double> &weights) {
    const auto size = cbls.size();

    for (auto i = 0ul; i < size; ++i) {
      const auto cbl_center = center(cbls[i]);

      const auto closest_cluster_center =
          center(*closest_cluster(first, last, cbl_center));

      const auto cbl_to_cluster_dist2 =
          (cbl_center - closest_cluster_center).squaredNorm();

      weights[i] = cbl_to_cluster_dist2;
    }
  }

};  // class Kmeans

/*! \brief provides the quality of the k-means guess for a given clustering.
 *
 * The k-means objective function is the sum of the squares of the distance of
 * the cluster's center to the center of each element in the cluster.
 *
 * \param cs a vector of clusters which must match the requirements for clusters
 * in Kmeans
 *
 * \see Kmeans
 */
template <typename Cluster>
double kmeans_objective(std::vector<Cluster> const &cs) {
  double sum = 0.0;
  // sum over clusters
  for (auto const &cluster : cs) {
    const auto cluster_center = center(cluster);

    // sum over squared distance of elem to cluster centers
    for (auto const &elem : cluster) {
      const auto elem_center = center(elem);
      sum += (elem_center - cluster_center).squaredNorm();
    }
  }

  return sum;
}

/*! @} */

}  // namespace clustering
}  // namespace mpqc

#endif  // MPQC_CLUSTERING_KMEANS_H
