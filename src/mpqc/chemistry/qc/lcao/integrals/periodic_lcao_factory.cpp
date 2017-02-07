#include "mpqc/chemistry/qc/lcao/integrals/periodic_lcao_factory.h"

namespace mpqc {
namespace lcao {
namespace detail {

float take_real(TiledArray::TensorD &result_tile, const TiledArray::TensorZ &arg_tile) {
    result_tile = TA::TensorD(arg_tile.range());

    const auto p0 = result_tile.range().lobound()[0];
    const auto pn = result_tile.range().upbound()[0];
    const auto q0 = result_tile.range().lobound()[1];
    const auto qn = result_tile.range().upbound()[1];
    const auto r0 = result_tile.range().lobound()[2];
    const auto rn = result_tile.range().upbound()[2];
    const auto s0 = result_tile.range().lobound()[3];
    const auto sn = result_tile.range().upbound()[3];

    auto tile_idx = 0;
    float norm = 0.0;
    for (auto p = p0; p < pn; ++p) {
        for (auto q = q0; q < qn; ++q) {
            for (auto r = r0; r < rn; ++r) {
                for (auto s = s0; s < sn; ++s, ++tile_idx) {
                    auto result_pqrs = arg_tile[tile_idx].real();
                    norm += std::abs(result_pqrs) * std::abs(result_pqrs);
                    result_tile[tile_idx] = result_pqrs;
                }
            }
        }
    }
    return std::sqrt(norm);
}

}  // namespace detail
}  // namespace lcao
}  // namespace mpqc
