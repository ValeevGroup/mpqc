

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_TASK_INTEGRALS_COMMON_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_TASK_INTEGRALS_COMMON_H_

#include <memory>
#include <array>
#include <vector>
#include <utility>

#include <tiledarray.h>


#include "mpqc/chemistry/qc/basis/basis.h"
#include "mpqc/chemistry/qc/integrals/integral_engine_pool.h"

namespace mpqc {
namespace integrals {

using Shell = mpqc::basis::Basis::Shell;
using ShellVec = std::vector<Shell>;

template <typename E>
using ShrPool = std::shared_ptr<mpqc::integrals::EnginePool<E>>;

template <unsigned long N>
using Barray = std::array<basis::Basis, N>;

using Bvector = std::vector<basis::Basis>;

template <typename T>
using OrdTileVec = std::vector<std::pair<unsigned long, T>>;

const static Shell unit_shell = Shell::unit();

namespace detail {

template <unsigned long N>
using VecArray = std::array<ShellVec const *, N>;

template <unsigned long N>
using ShrBases = std::shared_ptr<Barray<N>>;

using ShrBvectors = std::shared_ptr<Bvector>;

using IdxVec = std::vector<std::size_t>;

template <unsigned long N>
using ShrShellVecArray = std::array<std::shared_ptr<const ShellVec>, N>;

template <typename Op>
using Ttype = decltype(std::declval<Op>()(std::declval<TA::TensorD>()));

// Create TA::TiledRange from bases
template <unsigned long N>
TA::TiledRange create_trange(Barray<N> const &basis_array) {

    std::vector<TA::TiledRange1> trange1s;
    trange1s.reserve(N);

    for (auto i = 0ul; i < N; ++i) {
        trange1s.emplace_back(basis_array[i].create_trange1());
    }

    return TA::TiledRange(trange1s.begin(), trange1s.end());
}

// create TA::TiledRange from Bvector
inline TA::TiledRange create_trange(Bvector const& basis_vector) {

    std::size_t N = basis_vector.size();

    std::vector<TA::TiledRange1> trange1s;
    trange1s.reserve(N);

    for (auto i = 0ul; i < N; ++i) {
        trange1s.emplace_back(basis_vector[i].create_trange1());
    }

    return TA::TiledRange(trange1s.begin(), trange1s.end());

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
    auto const &pmap = a.pmap();
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

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_TASK_INTEGRALS_COMMON_H_
