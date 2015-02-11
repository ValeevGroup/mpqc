//
// task_integrals.h
//
// Copyright (C) 2014 Drew Lewis
// Maintainer Drew Lewis
//
// Based on file task_integrals.hpp from mpqc
//

#ifndef TCC_INTEGRALS_DENSETASKINTEGRALS_H
#define TCC_INTEGRALS_DENSETASKINTEGRALS_H

#include "../include/tiledarray.h"
#include "../include/btas.h"
#include "../include/libint2.h"

#include "../basis/cluster_shells.h"

#include "integral_engine_pool.h"
#include "btas_tensor_pass_through.h"
#include "tile_engine.h"

#include "../basis/basis.h"

#include "../tensor/btas_shallow_copy_wrapper.h"


namespace tcc {
namespace integrals {
    namespace dense{

template <std::size_t N>
using CShellPtrs
    = std::array<std::shared_ptr<std::vector<basis::ClusterShells>>, N>;

template <typename It, typename TF, typename SharedEnginePool,
          std::size_t N>
inline void do_task(It it, TiledArray::TiledRange trange, TF func,
            SharedEnginePool engines,
             CShellPtrs<N> shell_ptrs) {

    const auto idx = trange.tiles().idx(it.ordinal());
    const auto range = trange.make_tile_range(it.ordinal());

    const auto btas_tensor = tensor::ShallowTensor<N>{
        std::move(range), TileEngine<double>{}(idx, engines, shell_ptrs)};

    *it = func(std::move(btas_tensor));
}

template <typename Array, typename SharedEnginePool, std::size_t N, typename TF>
inline void
compute_tiles(Array &array, SharedEnginePool engines,
              std::array<basis::Basis, N> const &bases, TF tf) {

    CShellPtrs<N> shell_ptrs;

    for (auto i = 0ul; i < N; ++i) {
        shell_ptrs[i] = std::make_shared<std::vector<basis::ClusterShells>>(
            bases[i].cluster_shells());
    }

    const auto end = array.end();
    for (auto it = array.begin(); it != end; ++it) {
        array.get_world().taskq.add(do_task<decltype(it), TF, SharedEnginePool, N>,
                        it, array.trange(), tf, engines, shell_ptrs);
    }
    array.get_world().gop.fence();
}

/// Create a trange from the input bases.
template <std::size_t N>
inline TiledArray::TiledRange
create_trange(std::array<basis::Basis, N> const &basis_array) {

    std::vector<TiledArray::TiledRange1> trange1_collector;
    trange1_collector.reserve(N);

    for (auto i = 0ul; i < N; ++i) {
        trange1_collector.push_back(basis_array[i].create_flattend_trange1());
    }

    return TiledArray::TiledRange(trange1_collector.begin(),
                                  trange1_collector.end());
}
} // namespace dense

/// Create a dense array of tile where the type is specified by the functor.
/// TF needs to overload () to take a TiledArray::Range and a
/// btas::Tensor<double, btas::RangeNd<CblasRowMajor, std::array<long, N>>>;
template <typename SharedEnginePool, std::size_t N,
          typename TF = BtasTensorPassThrough<N>>
TiledArray::Array<double, N, typename TF::TileType, TiledArray::DensePolicy>
DenseIntegrals(madness::World &world, SharedEnginePool engines,
                     std::array<basis::Basis, N> const &bases, TF tf = TF{}) {

    auto trange = dense::create_trange(bases);
    TiledArray::Array<double, N, typename TF::TileType> array(world, trange);
    dense::compute_tiles(array, std::move(engines), bases, tf);

    return array;
}

} // namespace integrals
} // namespace tcc


#endif // TCC_INTEGRALS_DENSETASKINTEGRALS_H
