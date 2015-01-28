#pragma once
#ifndef TCC_INTEGRALS_BTASTENSORPASSTHROUGH_H
#define TCC_INTEGRALS_BTASTENSORPASSTHROUGH_H

#include "../include/btas.h"
#include "../tensor/btas_shallow_copy_wrapper.h"

namespace tcc {
namespace integrals {

/// This class is literally just a pass through for the default out type
/// of the tile_engine class.
template <unsigned int N>
class BtasTensorPassThrough {
  public:
    using TileType = tensor::ShallowTensor<N>;

    TileType operator()(TileType tile) const { return tile; }
};

} // namespace integrals
} // namespace tcc


#endif /* end of include guard: TCC_INTEGRALS_BTASTENSORPASSTHROUGH_H */
