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

} // namespace detail

} // namespace integrals
} // namespace mpqc

#endif //  MPQC_INTEGRALS_TASKINTEGRALSCOMMON_H
