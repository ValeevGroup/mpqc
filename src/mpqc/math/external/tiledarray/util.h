#ifndef MPQC4_SRC_MPQC_MATH_EXTERNAL_TILEDARRAY_UTIL_H_
#define MPQC4_SRC_MPQC_MATH_EXTERNAL_TILEDARRAY_UTIL_H_

#include <tiledarray.h>

namespace mpqc {
namespace detail {

/*!
 * \brief This extends 1D tiled range by repeating it multiple times
 *
 * \param tr0 the original TiledRange1 object
 * \param size the number of times for repeating
 * \return the extended TiledRange1 object
 */
TA::TiledRange1 extend_trange1(TA::TiledRange1 const &tr0, int64_t size);

}  // namespace detail
}  // namespace mpqc

#endif // MPQC4_SRC_MPQC_MATH_EXTERNAL_TILEDARRAY_UTIL_H_
