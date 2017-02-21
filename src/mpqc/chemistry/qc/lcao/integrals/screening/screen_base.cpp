/*
 * Created by Drew Lewis on 02/21/2017
 */

#include "screen_base.h"

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

TA::Tensor<float> Screener::norm_estimate(
    std::vector<gaussian::Basis> const &bs_array) const {

  auto trange = gaussian::detail::create_trange(bs_array);
  return TA::Tensor<float>(trange.tiles_range(),
                           std::numeric_limits<float>::max());

}

TA::Tensor<float> Screener::norm_estimate(madness::World &world,
    std::vector<gaussian::Basis> const &bs_array) const {

  return norm_estimate(bs_array);

}

}  // namespace lcao
}  // namespace mpqc
