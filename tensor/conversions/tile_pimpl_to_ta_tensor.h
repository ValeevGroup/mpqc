#pragma once
#ifndef TCC_TENSOR_CONVERSIONS_TILEPIMPLTOTATENSOR_H
#define TCC_TENSOR_CONVERSIONS_TILEPIMPLTOTATENSOR_H

#include "../../common/typedefs.h"

#include "../../include/tiledarray.h"
#include "../tile_pimpl.h"

namespace tcc {
namespace tensor {


template <typename T>
TA::Tensor<T> tile_pimpl_to_ta_tensor(TilePimpl<T> const &tp) {
    auto range = tp.range();
    return TA::Tensor<T>(tp.range(), tp.tile().matrix().data());
}

} // namespace tensor
} // namespace tcc

#endif // TCC_TENSOR_CONVERSIONS_TILEPIMPLTOTATENSOR_H
