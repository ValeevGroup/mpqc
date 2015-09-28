#pragma once

#ifndef MPQC_INTEGRALS_TASKINTEGRALSCOMMON_H
#define MPQC_INTEGRALS_TASKINTEGRALSCOMMON_H

#include "task_integrals_helper.h"
#include "../common/typedefs.h"
#include "../include/tiledarray.h"
#include "../include/tbb.h"
#include "../basis/basis.h"

#include <memory>
#include <array>

namespace mpqc {
namespace integrals {

template <typename E>
using ShrPool = std::shared_ptr<Epool<E>>;

template <unsigned long N>
using Barray = std::array<tcc::basis::Basis, N>;

namespace detail {

// Mutex for writing tiles to some shared structure
static tbb::spin_mutex task_int_mutex;

template <unsigned long N>
using ShrBases = std::shared_ptr<Barray<N>>;

template <typename Op>
using Ttype = decltype(std::declval<Op>()(std::declval<TA::TensorD>()));

// Create TRange from bases
template <unsigned long N>
TRange create_trange(Barray<N> const &basis_array) {

    std::vector<TRange1> trange1s;
    trange1s.reserve(N);

    for (auto i = 0ul; i < N; ++i) {
        trange1s.emplace_back(basis_array[i].create_flattend_trange1());
    }

    return TRange(trange1s.begin(), trange1s.end());
}

} // namespace detail

} // namespace integrals
} // namespace mpqc

#endif //  MPQC_INTEGRALS_TASKINTEGRALSCOMMON_H
