#ifndef MPQC_ARRAY_TILE_HPP
#define MPQC_ARRAY_TILE_HPP

#include "mpqc/range.hpp"
#include "mpqc/array/forward.hpp"

namespace mpqc {
namespace detail {

    struct ArrayTile {
	std::vector<range> extents;
	int proc, local;
	Array<void> *object;
	std::vector<range> subset(const std::vector<range> &ranges) const {
	    std::vector<range> s;
	    for (int i = 0; i < ranges.size(); ++i) {
		range r = range::intersection(ranges[i], this->extents[i]);
		if (!r.size()) return std::vector<range>();
		s.push_back(r);
	    }
	    return s;
	}
    };

} // namespace detail
} // namespace mpqc


#endif /* MPQC_ARRAY_TILE_HPP */
