#pragma once
#ifndef MPQC_CLUSTERING_COMMON_H
#define MPQC_CLUSTERING_COMMON_H

#include "../common/typedefs.h"

#include <vector>

namespace mpqc {
namespace clustering {

template <typename Iter>
Iter closest_cluster(Iter begin, Iter end, Vec3D const &cbls_center) {
    using Cluster = typename std::remove_reference<decltype(*begin)>::type;
    return std::min_element(begin, end,
                            [&](Cluster const &a, Cluster const &b) {
        return (center(a) - cbls_center).squaredNorm()
               < (center(b) - cbls_center).squaredNorm();
    });
}

template <typename Cluster>
double kmeans_objective(std::vector<Cluster> const &cs) {
    double sum = 0.0;
    // sum over clusters
    for(auto const &cluster : cs){
        const auto cluster_center = center(cluster);

        // sum over squared distance of elem to cluster centers
        for(auto const &elem : cluster){
            const auto elem_center = center(elem);
            sum += (elem_center - cluster_center).squaredNorm();
        }
    }

    return sum;
}

} // namespace clustering
} // namespace mpqc

#endif // MPQC_CLUSTERING_COMMON_H
