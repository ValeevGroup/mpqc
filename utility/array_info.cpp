#include "array_info.h"

#include "../tensor/decomposed_tensor.h"
#include "../tensor/mpqc_tile.h"

namespace mpqc {
namespace utility {

unsigned long
tile_clr_storage(Tile<DecompTensorD> const &tile) {
    auto size = 0ul;
    for (auto const &t : tile.tile().tensors()) {
        size += t.range().volume();
    }
    return size;
}

unsigned long tile_clr_storage(TA::TensorD const &) { return 0ul; }

} // namespace utility
} // namespace mpqc
