//
// task_integrals.h
//
// Copyright (C) 2014 Drew Lewis
// Maintainer Drew Lewis
//
// Based on file task_integrals.hpp from mpqc
//

#ifndef TCC_INTEGRALS_SPARSETASKINTEGRALS_H
#define TCC_INTEGRALS_SPARSETASKINTEGRALS_H

#include "../include/tiledarray.h"
#include "../include/btas.h"
#include "../include/libint2.h"

#include "../tensor/tile_algebra.h"

#include "../basis/cluster_shells.h"
#include "integral_engine_pool.h"
#include "btas_compute_functors.h"

#include "../basis/basis.h"

#include "../tensor/tile_pimpl_devel.h"

namespace tcc {
namespace integrals {

template <typename Pmap, typename SharedEnginePool, typename TileFunctor>
std::vector<std::array<double, 3>>
compute_tiles(madness::World &world, Pmap &p,
              TiledArray::TiledRange const &trange, SharedEnginePool engines,
              basis::Basis const &basis, TileFunctor tf) {

    std::vector<std::array<double, 3>> tiles(p->size());
    auto shared_tiles = std::make_shared<decltype(tiles)>(std::move(tiles));

    auto call_tf = [=](std::shared_ptr<decltype(tiles)> tile_vec_ptr,
                       TiledArray::Range range, TiledArray::Range::index idx,
                       basis::Basis const *basis_) {

        auto &vec = *tile_vec_ptr;
        auto const ordinal = trange.tiles().ordinal(idx);
        auto tile = tf(std::move(range), idx, basis_, engines);

        auto const &extent = tile.range().size();
        auto ij = extent[0] * extent[1];
        auto kl = extent[2] * extent[3];

        double dense_storage = double(ij * kl) * 8.0 / 1.0e9;
        double sparse_storage = dense_storage;
        double low_rank_storage = dense_storage;

        auto norm = tile.norm();
        if (norm < TiledArray::SparseShape<float>::threshold()) {
            sparse_storage = 0.0;
            low_rank_storage = 0.0;
        } else {
            Eigen::MatrixXd eig_tile = TiledArray::eigen_map(tile, ij, kl);
            Eigen::MatrixXd L, R;
            bool full_rank = algebra::Decompose_Matrix(eig_tile, L, R, 1e-8);

            if (!full_rank) {
                low_rank_storage = double(L.size() + R.size()) * 8.0 / 1.0e9;
                if (low_rank_storage > dense_storage) {
                    std::cout << "WTF MATE" << std::endl;
                    std::cout << "\tDense was " << dense_storage << " LR was "
                              << low_rank_storage << std::endl;
                    std::cout << "\tL.size() = " << L.size()
                              << " R.size() = " << R.size() << std::endl;
                    std::cout << "\tij = " << ij << " kl = " << kl
                              << " ij * kl = " << ij *kl << std::endl;
                }
            }
        }

        vec[ordinal] = {{low_rank_storage, sparse_storage, dense_storage}};
    };

    madness::Range<decltype(p->begin())> m_range{p->begin(), p->end()};

    world.taskq.for_each(m_range, [=](decltype(p->begin()) it) {
                             call_tf(shared_tiles, trange.make_tile_range(*it),
                                     trange.tiles().idx(*it), &basis);
                             return true;
                         }).get();

    return *shared_tiles;
}

template <std::size_t order>
TiledArray::TiledRange get_trange(basis::Basis const &basis) {
    std::vector<TiledArray::TiledRange1> trange1_collector;
    trange1_collector.reserve(order);

    const auto trange1 = basis.create_flattend_trange1();
    for (auto i = 0ul; i < order; ++i) {
        trange1_collector.push_back(trange1);
    }

    return TiledArray::TiledRange(trange1_collector.begin(),
                                  trange1_collector.end());
}

template <unsigned long order>
std::shared_ptr<TiledArray::Pmap>
create_pmap(madness::World &world, TiledArray::TiledRange const &trange) {
    return TiledArray::SparsePolicy::default_pmap(world,
                                                  trange.tiles().volume());
}

template <typename SharedEnginePool,
          typename TileFunctor = compute_functors::BtasTileFunctor<double>>
void Compute_Storage(madness::World &world, SharedEnginePool engines,
                     basis::Basis const &basis,
                     TileFunctor tf = TileFunctor{}) {

    constexpr auto tensor_order = integrals::pool_order<decltype(engines)>();
    using TileType = typename TileFunctor::TileType;

    auto trange = get_trange<tensor_order>(basis);

    auto pmap = create_pmap<tensor_order>(world, trange);

    auto storage_info
        = compute_tiles(world, pmap, trange, std::move(engines), basis, tf);

    std::array<double, 3> out = {{0.0, 0.0, 0.0}};
    for (auto const &group : storage_info) {
        out[0] += group[0];
        out[1] += group[1];
        out[2] += group[2];
    }
    world.gop.sum(&out[0], 3);

    if (world.rank() == 0) {
        std::cout << "Low rank storage was " << out[0] << "GB\n"
                  << "Sparse storage was " << out[1] << "GB\n"
                  << "Dense storage was " << out[2] << "GB" << std::endl;
    }
}

} // namespace integrals
} // namespace tcc


#endif // TCC_INTEGRALS_SPARSETASKINTEGRALS_H
