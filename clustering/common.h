#pragma once
#ifndef MPQC_CLUSTERING_COMMON_H
#define MPQC_CLUSTERING_COMMON_H

#include "../common/typedefs.h"

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

} // namespace clustering
} // namespace mpqc

#endif // MPQC_CLUSTERING_COMMON_H
