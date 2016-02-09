#include "array_info.h"

#include "../tensor/decomposed_tensor.h"
#include "../tensor/mpqc_tile.h"

namespace mpqc {
namespace utility {

double
tile_clr_storage(Tile<DecompTensorD> const &tile) {
    auto size = 0.0;
    for (auto const &t : tile.tile().tensors()) {
        size += t.range().volume();
    }
    return size;
}

double tile_clr_storage(TA::TensorD const &tile) { return 0.0; }

} // namespace utility
} // namespace mpqc
