#include "clustering_functions.h"
#include "common.h"
#include "molecule.h"

#include "atom_based_cluster.h"
#include "atom_based_cluster_concept.h"

#include "../clustering/kmeans.h"
#include "../clustering/common.h"

#include <cassert>
#include <numeric>
#include <iostream>
#include <random>


namespace mpqc {
namespace molecule {

namespace {
using ABCbl = AtomBasedClusterable;
using ABCbls = std::vector<ABCbl>;

ABCbls convert_to_clusterable(std::vector<AtomBasedCluster> const &clusters) {
    ABCbls clusterables;
    clusterables.reserve(clusters.size());

    for (auto const &cluster : clusters) {
        clusterables.emplace_back(cluster);
    }

    return clusterables;
}

ABCbls attach_hydrogens(ABCbls const &clusterables) {
    assert(clusterables.size() != 0);

    std::vector<AtomBasedCluster> clusters;
    ABCbls hydrogens;
    for (auto const &clusterable : clusterables) {
        if (clusterable.charge() != 1) {
            clusters.emplace_back(clusterable);
        } else {
            hydrogens.emplace_back(clusterable);
        }
    }

    // Check that not all atoms where hydrogens, if so just return
    // clusterables
    if (clusters.size() == 0) {
        return clusterables;
    }

    for (auto cluster : clusters) {
        update_center(cluster);
    }

    // Attach the hydrogens to the closest cluster
    if (!hydrogens.empty()) {
        for (auto &&hydrogen : hydrogens) {
            auto closest = clustering::closest_cluster(
                  clusters.begin(), clusters.end(), center(hydrogen));
            attach_clusterable(*closest, std::move(hydrogen));
        }

        for (auto cluster : clusters) {
            update_center(cluster);
        }
    }

    return convert_to_clusterable(clusters);
}

} // anon namespace 

Molecule
attach_hydrogens_and_kmeans(ABCbls const &clusterables, int64_t nclusters) {
    auto h_attached_clusterables = attach_hydrogens(clusterables);

    // TODO finish tomorrow add seeds
    clustering::Kmeans kmeans;
    return Molecule(convert_to_clusterable(kmeans.cluster<AtomBasedCluster>(
          h_attached_clusterables, nclusters)));
}

} // namespace molecule
} // namespace mpqc
