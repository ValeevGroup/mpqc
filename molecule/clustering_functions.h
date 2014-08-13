#ifndef CLUSTERING_FUNCTIONS_H
#define CLUSTERING_FUNCTIONS_H

#include "cluster_concept.h"
#include "cluster.h"
#include <random>

namespace clustering {
using output_t = std::vector<Cluster>;
using input_t = std::vector<Clusterable>;

/**
* @brief sum_cluster_distances sums the distance of each Clusterable to its
* cluster and then reduces them.
* @param clusters the clusters over which the distances will be reduced.
* @return the reduced distances of the atoms for each cluster.
*/
double sum_cluster_distances(const std::vector<Cluster> &clusters);

/**
 * @brief The kmeans class is a functor which kmeans clusters clusterables.
 */
class kmeans {
public:
  /**
   * @brief kmeans constructor provides the seed to generate the random inital
   * guess.
   * @param seed is an unsigned type which is used to inialize the random number
   * generator
   * used for the inital cluster guess.
   */
  kmeans(unsigned long seed);

  /**
   * @brief operator () computes clusters over the clusterables
   * @param clusterables the elements to be clustered.
   * @param nclusters the number of clusters desired.
   * @return a vector of clusters.
   */
  output_t operator()(input_t clusterables, unsigned long nclusters);

private:
  /**
   * @brief initialize_clusters takes the clusterables which need to be
   * clustered
   * @param clusterables const reference to objects which need to be clustered.
   */
  void initialize_clusters(const input_t &clusterables);

  /**
   * @brief attach_clusterables loops over clusterables and attaches them to
   * the nearest cluster.
   * @param cs a ref to a const vector of clusterables
   */
  void attach_clusterables(const std::vector<Clusterable> &cs);

  /**
   * @brief update_clusters calculates new center based on old clusters.
   * @param clusterables is a ref to a const vector of clustersables
   * @return the previous iteration's centers.
   */
  std::vector<Cluster::position_t>
  update_clusters(const std::vector<Clusterable> &clusterables);

  /**
   * @brief kmeans_converged determines whether or not the kmeans calcualation
   * has converged.
   * @param old_centers is the previous iterations centers.
   */
  bool kmeans_converged(const std::vector<Cluster::position_t> &old_centers);

  /**
   * @brief cluster performs the kmeans iterations to generate clusters.
   * @param clusterables a const reference to a vector of clusterabls.
   * @return the clusters.
   */
  output_t cluster(const std::vector<Clusterable> &clusterables);


  std::vector<Cluster>::iterator
  closest_cluster(const std::vector<Cluster>::iterator begin,
                  const std::vector<Cluster>::iterator end, const position_t center);

  unsigned long seed_;
  std::vector<Cluster> clusters_;
};

} // namespace clustering

#endif // CLUSTERING_FUNCTIONS_H
