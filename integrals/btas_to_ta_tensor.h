#pragma once
#ifndef TCC_INTEGRALS_BTASCOMPUTEFUNCTORS_H
#define TCC_INTEGRALS_BTASCOMPUTEFUNCTORS_H

#include "../include/btas.h"
#include "../include/tiledarray.h"

#include "tile_engine.h"

namespace tcc {
namespace integrals {
namespace compute_functors {

struct BtasToTaTensor {
    using Tensor = TiledArray::Tensor<double>;
    using TileType = Tensor;

    template <std::size_t N>
    Tensor operator()(TiledArray::Range const &range,
                      TileEngine<double>::TileType<N> const &btas_tensor) {
        Tensor tensor{range};

        const auto b_data = btas_tensor.data();
        const auto size = btas_tensor.size();
        std::copy(b_data, b_data + size, tensor.data());

        return tensor;
    }
}; // BtasToTaTensor


} // namespace compute_functors
} // namespace integrals
} // namespace tcc

#endif /* end of include guard: TCC_INTEGRALS_BTASCOMPUTEFUNCTORS_H */
