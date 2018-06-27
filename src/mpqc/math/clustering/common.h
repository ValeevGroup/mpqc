
#ifndef MPQC4_SRC_MPQC_MATH_CLUSTERING_COMMON_H_
#define MPQC4_SRC_MPQC_MATH_CLUSTERING_COMMON_H_

#include <vector>

namespace mpqc {
namespace math {
namespace clustering {

/*! \ingroup MathClustering @{ */

/*! \brief returns the closest element in the input range to the center provided
 *
 *  \param begin an iterator to the first cluster in the list to search over.
 *  \param end an iterator to one past the last cluster to search over.
 *  \param target_center a Vector3d, the center, for which the closest cluster
 *      will be found
 *
 *  The cluster type only requires that a non-intrusive center function that
 *  returns a Vec3d be defined.
 */
template <typename Iter>
Iter closest_cluster(Iter begin, Iter end, Vector3d const &target_center) {
  using Cluster = typename std::remove_reference<decltype(*begin)>::type;
  return std::min_element(begin, end, [&](Cluster const &a, Cluster const &b) {
    return (center(a) - target_center).squaredNorm() <
           (center(b) - target_center).squaredNorm();
  });
}

/*! @} */

}  // namespace clustering
}  // namespace math
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_MATH_CLUSTERING_COMMON_H_
