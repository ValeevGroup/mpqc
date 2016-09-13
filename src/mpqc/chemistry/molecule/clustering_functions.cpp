#include <mpqc/chemistry/molecule/clustering_functions.h>
#include <mpqc/chemistry/molecule/common.h>
#include <mpqc/chemistry/molecule/molecule.h>

#include "../../../../clustering/kmeans.h"


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

        for (auto &cluster : clusters) {
            update_center(cluster);
        }
    }

    return convert_to_clusterable(clusters);
}

} // anon namespace

Molecule
attach_hydrogens_and_kmeans(ABCbls const &clusterables, int64_t nclusters) {
    auto h_attached_clusterables = attach_hydrogens(clusterables);
    return kmeans(h_attached_clusterables, nclusters);
}

Molecule kmeans(ABCbls const &clusterables, int64_t nclusters) {
    if (clusterables.size() < std::size_t(nclusters)) {
        std::cout << "\nWarning!! User asked for more clusters than there were clusterables! Use "<<  clusterables.size() << " clusters!\n" << std::endl;
        nclusters = clusterables.size();
    }

    auto objective_min = std::numeric_limits<double>::max();
    int64_t init_seed = 1000;
    int64_t best_seed = 1000;
    for (auto i = 0; i < 50; ++i) {
        clustering::Kmeans kmeans(init_seed);
        auto clusters
              = kmeans.cluster<AtomBasedCluster>(clusterables, nclusters);
        auto value = clustering::kmeans_objective(clusters);
        if (value < objective_min) {
            best_seed = init_seed;
            objective_min = value;
        }
        init_seed += 1000;
    }

    return Molecule(convert_to_clusterable(
          clustering::Kmeans(best_seed)
                .cluster<AtomBasedCluster>(clusterables, nclusters)));
}

} // namespace molecule
} // namespace mpqc
