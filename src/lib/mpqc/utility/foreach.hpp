#ifndef MPQC_FOREACH_HPP
#define MPQC_FOREACH_HPP

#include <boost/foreach.hpp>

/// foreach loop construct, short-hand for BOOST_FOREACH
/// @use <c>foreach(const auto &r : range(0,1)) {}</c>
/// @ingroup Utility
#define foreach BOOST_FOREACH

#endif /* MPQC_FOREACH_HPP */
