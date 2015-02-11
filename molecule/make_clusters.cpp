#include "cluster.h"
#include "molecule.h"
#include "make_clusters.h"

namespace tcc {
namespace molecule {
std::vector<std::shared_ptr<Cluster>> attach_hydrogens_kmeans(Molecule const &m, unsigned long nclusters){
    std::vector<std::shared_ptr<Cluster>> clusters;
    clusters.reserve(nclusters);

    for(auto &&cluster : m.attach_H_and_kmeans(nclusters)){
        clusters.push_back(std::make_shared<Cluster>(std::move(cluster)));
    }
    
    return clusters;
}

} // namespace molecule
} // namespace tcc
