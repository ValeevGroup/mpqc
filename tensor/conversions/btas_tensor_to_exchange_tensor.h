#pragma once

#include "../../common/namespaces.h"
#include "../../common/typedefs.h"
#include "../../include/btas.h"
#include "../../include/tiledarray.h"

#include "../btas_shallow_copy_wrapper.h"
#include "../decomp_tensor.h"
#include "../coulomb_tensor.h"
#include "../decomp_algebra.h"

namespace tcc {
namespace tensor {
namespace conversions {

struct BtasToExchange {
    using TileType = ExchangeTensor;

    double cut;
    BtasToExchange(double c) : cut(c) {}

    ExchangeTensor operator()(ShallowTensor<4> const &bt) {
        if (bt.norm() < TA::SparseShape<float>::threshold()) {
            return ExchangeTensor{bt.range(), TATensor{}, cut};
        } else {
            // create btas view
            auto k_range = btas::permute(bt.tensor().range(), {0, 3, 2, 1});
            remove_ref_const_t<decltype(bt.tensor())> view
                = btas::make_view(k_range, bt.tensor().storage());

            auto const &extent = k_range.extent();
            TARange range{extent[0], extent[1], extent[2], extent[3]};
            TATensor t(range);

            std::copy(view.begin(), view.end(), t.data());

            auto tuple = algebra::exchange_decomposition(t, cut);
            return ExchangeTensor(tuple, cut);
        }
    }
};

struct BtasToCoulomb {
    using TileType = CoulombTensor;

    double cut;
    BtasToCoulomb(double c) : cut(c) {}

    CoulombTensor operator()(ShallowTensor<4> const &bt) {
        if (bt.norm() < TA::SparseShape<float>::threshold()) {
            return CoulombTensor{bt.range(), TATensor{}, cut};
        } else {
            TARange range{bt.range()};
            TATensor t(range);

            std::copy(bt.tensor().begin(), bt.tensor().end(), t.data());

            auto tuple = algebra::coulomb_decomposition(t, cut);
            return CoulombTensor(tuple, cut);
        }
    }
};

} // namespace conversions
} // namespace tensor
} // namespace tcc
