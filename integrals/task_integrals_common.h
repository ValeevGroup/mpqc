#pragma once

#ifndef MPQC_INTEGRALS_TASKINTEGRALSCOMMON_H
#define MPQC_INTEGRALS_TASKINTEGRALSCOMMON_H

#include "../common/typedefs.h"
#include "../basis/basis.h"

#include "../include/tiledarray.h"

#include <memory>
#include <array>
#include <vector>
#include <utility>

namespace mpqc {
namespace integrals {

template <typename E>
using ShrPool = std::shared_ptr<Epool<E>>;

template <unsigned long N>
using Barray = std::array<basis::Basis, N>;

template <typename T>
using OrdTileVec = std::vector<std::pair<unsigned long, T>>;

const static Shell unit_shell = Shell::unit();

struct TensorPassThrough {
    TA::TensorD operator()(TA::TensorD &&ten) const {
        return std::move(ten);
    }
};

namespace detail {

template <unsigned long N>
using VecArray = std::array<ShellVec const *, N>;

template <unsigned long N>
using ShrBases = std::shared_ptr<Barray<N>>;

using IdxVec = std::vector<std::size_t>;

template <unsigned long N>
using ShrShellVecArray = std::array<std::shared_ptr<const ShellVec>, N>;

template <typename Op>
using Ttype = decltype(std::declval<Op>()(std::declval<TA::TensorD>()));

// Create TRange from bases
template <unsigned long N>
TRange create_trange(Barray<N> const &basis_array) {

    std::vector<TRange1> trange1s;
    trange1s.reserve(N);

    for (auto i = 0ul; i < N; ++i) {
        trange1s.emplace_back(basis_array[i].create_trange1());
    }

    return TRange(trange1s.begin(), trange1s.end());
}

template<unsigned long N>
ShrShellVecArray<N> get_shells(IdxVec const &idx, ShrBases<N> const& bases){

    ShrShellVecArray<N> shell_vecs;
    for(auto i = 0ul; i < N; ++i){
        auto const &basis = (*bases)[i];
        auto *shell_vec = &(basis.cluster_shells()[idx[i]]);
        shell_vecs[i] = std::shared_ptr<const ShellVec>(bases, shell_vec);
    }

    return shell_vecs;
}

template<typename Tile, typename Array>
void set_array(std::vector<std::pair<unsigned long, Tile>> &tiles, Array &a){
    auto const &pmap = a.get_pmap();
    for(auto ord : *pmap){
        if(a.is_local(ord) && !a.is_zero(ord)){
            auto &tile = tiles[ord];
            assert(!tile.second.empty());
            a.set(ord, std::move(tile.second));
        }
    }
}

} // namespace detail

} // namespace integrals
} // namespace mpqc

#endif //  MPQC_INTEGRALS_TASKINTEGRALSCOMMON_H
