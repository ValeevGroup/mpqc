#pragma once
#ifndef MPQC_CLUSTERING_COMMON_H
#define MPQC_CLUSTERING_COMMON_H

#include "../common/typedefs.h"

#include <vector>

namespace mpqc {
namespace clustering {

/*! \ingroup Clustering @{ */


/*! \brief returns the closest element in the input range to the center provided
 *
 *  \param begin an iterator to the first cluster in the list to search over.
 *  \param end an iterator to one past the last cluster to search over.
 *  \param target_center a Vec3D, the center, for which the closest cluster
 *      will be found 
 *
 *  The cluster type only requires that a non-intrusive center function that
 *  returns a Vec3d be defined.
 */
template <typename Iter>
Iter closest_cluster(Iter begin, Iter end, Vec3D const &target_center) {
    using Cluster = typename std::remove_reference<decltype(*begin)>::type;
    return std::min_element(begin, end,
                            [&](Cluster const &a, Cluster const &b) {
        return (center(a) - target_center).squaredNorm()
               < (center(b) - target_center).squaredNorm();
    });
}

/*! @} */


} // namespace clustering
} // namespace mpqc

#endif // MPQC_CLUSTERING_COMMON_H
