#include "mpqc/math/external/tiledarray/util.h"

namespace mpqc {
namespace detail {

TA::TiledRange1 extend_trange1(TA::TiledRange1 const &tr0, int64_t size) {
  auto blocking = std::vector<int64_t>{0};
  for (decltype(size) idx = 0; idx < size; ++idx) {
    for (decltype(tr0.tile_extent()) u = 0; u < tr0.tile_extent(); ++u) {
      auto next = blocking.back() + tr0.tile(u).second - tr0.tile(u).first;
      blocking.emplace_back(next);
    }
  }
  TA::TiledRange1 tr1(blocking.begin(), blocking.end());
  return tr1;
}

}  // namespace detail
}  // namespace mpqc

