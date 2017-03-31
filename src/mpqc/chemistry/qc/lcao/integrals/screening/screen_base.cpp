/*
 * Created by Drew Lewis on 02/21/2017
 */

#include "./screen_base.h"

#include "mpqc/util/misc/exception.h"

namespace mpqc {
namespace lcao {

bool Screener::skip(int64_t) { return false; }
bool Screener::skip(int64_t) const { return false; }

bool Screener::skip(int64_t, int64_t) { return false; }
bool Screener::skip(int64_t, int64_t) const { return false; }

bool Screener::skip(int64_t, int64_t, int64_t) { return false; }
bool Screener::skip(int64_t, int64_t, int64_t) const { return false; }

bool Screener::skip(int64_t, int64_t, int64_t, int64_t) { return false; }
bool Screener::skip(int64_t, int64_t, int64_t, int64_t) const { return false; }

bool Screener::skip(int64_t, int64_t, int64_t, int64_t, double) { return false; }
bool Screener::skip(int64_t, int64_t, int64_t, int64_t, double) const { return false; }

TA::Tensor<float> Screener::norm_estimate(
    madness::World &world, std::vector<gaussian::Basis> const &bs_array,
    TA::Pmap const &pmap, bool replicate) const {
  const auto trange = gaussian::detail::create_trange(bs_array);
  const auto ndim = bs_array.size();
  if (ndim == 2) {
    assert(false && "not implemented");
  } else if (ndim == 3) {
    assert(false && "not implemented");
  } else if (ndim == 4) {
    assert(false && "not implemented");
  } else {
    throw ProgrammingError(
        "Screener::norm_estimate only supports 2, 3, and 4 dimensions",
        __FILE__, __LINE__);
  }

  return TA::Tensor<float>(trange.tiles_range(),
                           std::numeric_limits<float>::max());
}

}  // namespace lcao
}  // namespace mpqc
