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

#include "../basis/cluster_shells.h"

#include "integral_engine_pool.h"
#include "btas_tensor_pass_through.h"
#include "tile_engine.h"

#include "../basis/basis.h"

#include "../tensor/btas_shallow_copy_wrapper.h"

#include "../common/namespaces.h"


namespace tcc {
namespace integrals {
namespace sparse {

template <std::size_t N>
using CShellPtrs
      = std::array<std::shared_ptr<std::vector<basis::ClusterShells>>, N>;

template <typename TileType, typename TF, typename SharedEnginePool,
          std::size_t N>
void do_task(std::vector<std::pair<std::size_t, TileType>> *tile_vec,
             std::size_t ord, std::size_t tile_ord, TF func,
             TiledArray::TiledRange trange, SharedEnginePool engines,
             CShellPtrs<N> shell_ptrs) {
    auto &vec = *tile_vec;
    const auto idx = trange.tiles().idx(tile_ord);
    auto range = trange.make_tile_range(tile_ord);

    /*! \todo Fix this to return the correct type and avoid copies */
    const auto btas_tensor = tensor::ShallowTensor<N>{
          std::move(range), TileEngine<double>{}(idx, engines, shell_ptrs)};

    auto ta_tensor = func(btas_tensor);
    auto norm = ta_tensor.norm();
    assert(norm >= 0);
    assert(norm == norm);

    // Save the tensor with it's ordinal information for later.
    vec[ord] = std::make_pair(tile_ord, func(std::move(btas_tensor)));
}

template <typename TileType, typename SharedEnginePool, std::size_t N>
void do_multipole_task(
      std::vector<std::vector<std::pair<std::size_t, TileType>>> *tile_vec,
      std::size_t ord, std::size_t tile_ord, TiledArray::TiledRange trange,
      SharedEnginePool engines, CShellPtrs<N> shell_ptrs,
      unsigned long number_of_arrays) {

    auto &vec = *tile_vec;
    const auto idx = trange.tiles().idx(tile_ord);
    auto range = trange.make_tile_range(tile_ord);
    const auto extent = range.extent();

    std::vector<TileType*> tile_ptrs(number_of_arrays); 
    for (auto i = 0ul; i < number_of_arrays; ++i) {
        vec[i][ord] = std::make_pair(tile_ord, trange.make_tile_range(tile_ord));
        tile_ptrs[i] = &(vec[i][ord].second);
    }

    // Only work for matrices for now
    auto &shells0 = *shell_ptrs[0];
    auto &shells1 = *shell_ptrs[1];

    auto &engine = engines->local();
    using Mat = Eig::Matrix<double, Eig::Dynamic, Eig::Dynamic, Eig::RowMajor>;

    auto bf0 = 0;
    for (auto &sh0 : shells0[idx[0]].flattened_shells()) {
        auto n0 = sh0.size();
        auto bf1 = 0;
        for (auto &sh1 : shells1[idx[1]].flattened_shells()) {
            auto n1 = sh1.size();
            auto n01 = n0 * n1;

            const auto *buf = engine.compute(sh0, sh1);

            for(auto i = 0ul; i < number_of_arrays; ++i){
                auto buf_map = Eig::Map<const Mat>(buf + (i * n01), n0, n1);
                auto tmap = TA::eigen_map(*tile_ptrs[i], extent[0], extent[1]);

                tmap.block(bf0, bf1, n0, n1) = buf_map;
            }

            bf1 += n1;
        }
        bf0 += n0;
    }
}


template <typename Pmap, typename SharedEnginePool, std::size_t N, typename TF>
std::vector<std::pair<std::size_t, typename TF::TileType>>
compute_tiles(madness::World &world, Pmap const &p,
              TiledArray::TiledRange const &trange, SharedEnginePool engines,
              std::array<basis::Basis, N> const &bases, TF tf) {

    CShellPtrs<N> shell_ptrs;

    for (auto i = 0ul; i < N; ++i) {
        shell_ptrs[i] = std::make_shared<std::vector<basis::ClusterShells>>(
              bases[i].cluster_shells());
    }

    const auto begin = p->begin();
    const auto end = p->end();
    const auto ntiles = std::distance(begin, end);
    std::vector<std::pair<std::size_t, typename TF::TileType>> tiles(ntiles);

    auto ord = 0ul;
    for (auto it = begin; it != end; ++it, ++ord) {
        const auto idx = trange.tiles().idx(*it);
        world.taskq.add(do_task<typename TF::TileType, TF, SharedEnginePool, N>,
                        &tiles, ord, *it, tf, trange, engines, shell_ptrs);
    }
    world.gop.fence();

    return tiles;
}

template <typename Pmap, typename SharedEnginePool, std::size_t N, typename TF>
std::vector<std::vector<std::pair<std::size_t, typename TF::TileType>>>
compute_tiles(madness::World &world, Pmap const &p,
              TiledArray::TiledRange const &trange, SharedEnginePool engines,
              std::array<basis::Basis, N> const &bases, 
              unsigned long number_of_arrays, TF const &tf) {

    CShellPtrs<N> shell_ptrs;

    for (auto i = 0ul; i < N; ++i) {
        shell_ptrs[i] = std::make_shared<std::vector<basis::ClusterShells>>(
              bases[i].cluster_shells());
    }

    const auto begin = p->begin();
    const auto end = p->end();
    const auto ntiles = std::distance(begin, end);

    using tile_vec = std::vector<std::pair<std::size_t, typename TF::TileType>>;
    std::vector<tile_vec> tiles;
    tiles.reserve(number_of_arrays);
    for (auto i = 0ul; i < number_of_arrays; ++i) {
        tiles.push_back(tile_vec(ntiles));
    }

    auto ord = 0ul;
    for (auto it = begin; it != end; ++it, ++ord) {
        const auto idx = trange.tiles().idx(*it);
        world.taskq.add(
              do_multipole_task<typename TF::TileType, SharedEnginePool, N>,
              &tiles, ord, *it, trange, engines, shell_ptrs, number_of_arrays);
    }
    world.gop.fence();

    return tiles;
}

/// Create a trange from the input bases.
template <std::size_t N>
TiledArray::TiledRange
create_trange(std::array<basis::Basis, N> const &basis_array) {

    std::vector<TiledArray::TiledRange1> trange1_collector;
    trange1_collector.reserve(N);

    for (auto i = 0ul; i < N; ++i) {
        trange1_collector.push_back(basis_array[i].create_flattend_trange1());
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

/// Creates a sparse array given a vector of pairs.  The pairs are ordinals with
/// the corresponding tile.
template <std::size_t N, typename TileType>
TiledArray::Array<double, N, TileType, TiledArray::SparsePolicy> create_array(
      madness::World &world, TiledArray::TiledRange const &trange,
      std::shared_ptr<TiledArray::Pmap> const &pmap,
      std::vector<std::pair<std::size_t, TileType>> const &ord_tile_pairs) {

    TiledArray::Tensor<float> tile_norms(trange.tiles(), 0.0);
    for (auto const &pair : ord_tile_pairs) {
        if (!pair.second.empty()) {
            tile_norms[pair.first] = pair.second.norm();
        }
    }

    TiledArray::SparseShape<float> shape(world, tile_norms, trange);

    TiledArray::Array<double, N, TileType, TiledArray::SparsePolicy> array(
          world, trange, shape, pmap);

    for (auto const &pair : ord_tile_pairs) {
        if (!array.is_zero(pair.first)) {
            array.set(pair.first, pair.second);
        }
    }

    return array;
}

} // namespace sparse

/// Create a sparse array of tile where the type is specified by the functor.
/// TF needs to overload () to take a TiledArray::Range and a
/// btas::Tensor<double, btas::RangeNd<CblasRowMajor, std::array<long, N>>>;
template <typename SharedEnginePool, std::size_t N,
          typename TF = BtasTensorPassThrough<N>>
TiledArray::Array<double, N, typename TF::TileType, TiledArray::SparsePolicy>
BlockSparseIntegrals(madness::World &world, SharedEnginePool engines,
                     std::array<basis::Basis, N> const &bases, TF tf = TF{}) {

    auto trange = sparse::create_trange(bases);
    auto pmap = sparse::create_pmap<N>(world, trange);

    // Compute the tiles for the array and save them in a vector
    auto computed_tiles = sparse::compute_tiles(world, pmap, trange,
                                                std::move(engines), bases, tf);

    return sparse::create_array<N>(world, trange, pmap, computed_tiles);
}

template <std::size_t N, typename TF = BtasTensorPassThrough<N>>
std::vector<TiledArray::Array<double, N, typename TF::TileType,
                              TiledArray::SparsePolicy>>
BlockSparseIntegrals(
      madness::World &world,
      std::shared_ptr<EnginePool<libint2::OneBodyEngine>> &engines,
      std::array<basis::Basis, N> const &bases, TF tf,
      unsigned long number_of_arrays) {

    auto trange = sparse::create_trange(bases);
    auto pmap = sparse::create_pmap<N>(world, trange);

    // Compute the tiles for the array and save them in a vector
    auto computed_tiles
          = sparse::compute_tiles(world, pmap, trange, std::move(engines),
                                  bases, number_of_arrays, tf);

    std::vector<TA::Array<double, N, typename TF::TileType, TA::SparsePolicy>>
          oxyz;
    for (auto i = 0ul; i < number_of_arrays; ++i) {
        oxyz.push_back(
              sparse::create_array<N>(world, trange, pmap, computed_tiles[i]));
    }
    return oxyz;
}


} // namespace integrals
} // namespace tcc


#endif // TCC_INTEGRALS_SPARSETASKINTEGRALS_H
