#ifndef MPQC_UTILITY_FOREACH_HPP
#define MPQC_UTILITY_FOREACH_HPP

#include <boost/foreach.hpp>

/// @ingroup CoreUtility
/// @brief foreach loop construct, short-hand for BOOST_FOREACH.
/// Example:
/// @code
/// foreach (const auto &r : range(0,1)) {}
/// @endcode
#define foreach BOOST_FOREACH

#endif /* MPQC_UTILITY_FOREACH_HPP */
