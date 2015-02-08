#pragma once
#ifndef TCC_INTEGRALS_LAZYINTEGRALS_H
#define TCC_INTEGRALS_LAZYINTEGRALS_H

#include "../include/tiledarray.h"
#include "../tensor/tile_pimpl_devel.h"
#include "integral_engine_pool.h"
#include "../basis/cluster_shells.h"
#include "../tensor/tile_variant_devel.h"
#include "../include/btas.h"
#include "btas_compute_functors.h"

#include <memory>

namespace tcc {
namespace integrals {

template <typename SharedEnginePool>
class LazyTile {
  public:
    using eval_type = tensor::TilePimplDevel<double, TiledArray::Range>;
    using value_type = eval_type;

    LazyTile() = default;
    LazyTile(LazyTile const &other) = default;
    LazyTile &operator=(LazyTile const &other) = default;
    LazyTile(LazyTile &&other) = default;
    LazyTile &operator=(LazyTile &&other) = default;

    LazyTile(TiledArray::Range range,
             SharedEnginePool engines_ptr,
             std::vector<basis::ClusterShells> clustered_shells)
        : range_{std::move(range)},
          engines_ptr_{std::move(engines_ptr)},
          clustered_shells_{std::move(clustered_shells)} {}

    operator eval_type() const {
        btas::Tensor<double> btas_tensor = compute_functors::detail::
            integral_wrapper<double, pool_order<decltype(engines_ptr_)>()>{}(
                engines_ptr_->local(), clustered_shells_);
        return eval_type(
            range_, tensor::TileVariantDevel<double>{std::move(btas_tensor)});
    }


    template <typename Archive>
    void serialize(const Archive &) {
        assert(false);
    }

  private:
    TiledArray::Range range_;
    SharedEnginePool engines_ptr_;
    std::vector<basis::ClusterShells> clustered_shells_;
}; // LazyTile

template <typename SharedEnginePool>
TiledArray::Array<double, pool_order<SharedEnginePool>(),
                  LazyTile<SharedEnginePool>>
create_lazy_array(madness::World &world, basis::Basis const &basis,
                  SharedEnginePool engines_ptr) {
    constexpr auto order = pool_order<SharedEnginePool>();

    std::vector<TiledArray::TiledRange1> trange1s;
    trange1s.reserve(order);

    const auto trange1 = basis.create_flattend_trange1();
    for (auto i = 0ul; i < order; ++i) {
        trange1s.push_back(trange1);
    }

    TiledArray::TiledRange trange(trange1s.begin(), trange1s.end());

    TiledArray::Array<double, order, LazyTile<SharedEnginePool>> array(
        world, trange);

    for (auto it = array.begin(); it != array.end(); ++it) {
        auto idx = it.index();
        std::vector<basis::ClusterShells> cluster_shells;
        cluster_shells.reserve(idx.size());
        for (auto i : idx) {
            cluster_shells.push_back(basis.cluster_shells()[i]);
        }

        *it = LazyTile<SharedEnginePool>{it.make_range(), engines_ptr,
                                           std::move(cluster_shells)};
    }

    return array;
}


} // integrals
} // tcc


#endif /* end of include guard: TCC_INTEGRALS_LAZYINTEGRALS_H */
