#pragma once
#ifndef TCC_UTLITY_TAHELPERS_H
#define TCC_UTLITY_TAHELPERS_H

#include "../common/namespaces.h"
#include "../common/typedefs.h"
#include "../include/tiledarray.h"

#include <string>

namespace mpqc {
namespace utility {

namespace detail {
static constexpr char letters[] = "abcdefghij";

// Helper function to make strings of indices, not particularly efficient, but 
// oh well, if string ever gets constexpr constructor comeback and fix.
template <unsigned int N>
std::string make_indices() {
    return std::string{letters[N-1], ','} + make_indices<N - 1>();
}

// Specialization to end the recursion. 
template <>
std::string make_indices<1u>() {
    return std::string{letters[0]};
}
}

// Function which will compute the Frobenius norm of the difference of
// two tensors which share a tile type and a policy.
template <typename T, unsigned int DIM, typename TileType, typename Policy>
double array_fnorm_diff(TA::Array<T, DIM, TileType, Policy> const &A,
                        TA::Array<T, DIM, TileType, Policy> const &B) {
    static_assert(DIM >= 1, "The dimension must be at least 1");
    static_assert(DIM < 10, "Dimensions above 10 are not currently supported");

    // Make output array and subtract A and B;
    remove_ref_const_t<decltype(A)> out;
    const std::string indices = detail::make_indices<DIM>();
    out(indices) = A(indices) - B(indices);
    return out(indices).norm(out.get_world()).get();
}

} // namespace utility
} // namespace mpqc

#endif // TCC_UTLITY_TAHELPERS_H
