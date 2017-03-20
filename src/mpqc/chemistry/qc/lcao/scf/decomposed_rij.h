#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_DECOMPOSED_RIJ_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_DECOMPOSED_RIJ_H_

#include "mpqc/chemistry/qc/lcao/expression/orbital_index.h"
#include "mpqc/chemistry/qc/lcao/basis/basis.h"
#include "mpqc/util/misc/exception.h"
#include "mpqc/math/tensor/clr/array_to_eigen.h"
#include "mpqc/util/misc/exenv.h"

namespace mpqc {
namespace lcao {
namespace gaussian {

namespace detail {

template <typename T>
using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

/// create a unit basis that only has a unit shell
std::shared_ptr<Basis> create_unit_basis();

/// take \c row from a 2-d array and form a 1-d array
template<typename Tile, typename Policy, typename TArray = TA::DistArray<Tile, Policy>>
TArray take_row_from_2D_array(TA::DistArray<Tile, Policy> & matrix, int64_t row_idx) {
    using element_type = typename TArray::element_type;

    auto tr0 = matrix.trange().data()[0];
    auto tr1 = matrix.trange().data()[1];

    auto nrows = tr0.extent();

    if (row_idx < 0 || row_idx >= nrows)
        throw ProgrammingError("Requested row index is out of range", __FILE__, __LINE__);

    auto mat_eig = array_ops::array_to_eigen(matrix);
    Vector<element_type> vec_eig = mat_eig.row(row_idx);

    auto& world = matrix.world();
    TA::TiledRange trange{tr1};
    TA::TensorF norms(trange.tiles_range(), trange.elements_range().volume() / trange.tiles_range().volume());
    TA::SparseShape<float> shape(world, norms, trange);

    TArray result(world, trange, shape);

    auto vec_to_tile = [](TA::Range &range, const Vector<element_type>& arg_vec){
        auto tile = Tile(range);
        const auto p0 = tile.range().lobound()[0];
        const auto pn = tile.range().upbound()[0];

        auto tile_ord = 0;
        for (auto p = p0; p < pn; ++p, ++tile_ord) {
            const auto element = arg_vec[p];
            tile[tile_ord] = element;
        }

        return tile;
    };

    auto const &pmap = result.pmap();
    for (auto it = pmap->begin(); it != pmap->end(); ++it) {
        if (!result.is_zero(*it)) {
            auto range = trange.make_tile_range(*it);
            auto tile = world.taskq.add(vec_to_tile, range, vec_eig);
            result.set(*it, tile);
        }
    }
    world.gop.fence();
    result.truncate();
    return result;
}

/// take \c col from a 2-d array and form a 1-d array
template<typename Tile, typename Policy, typename TArray = TA::DistArray<Tile, Policy>>
TArray take_col_from_2D_array(TA::DistArray<Tile, Policy> & matrix, int64_t col_idx) {
    using element_type = typename TArray::element_type;

    auto tr0 = matrix.trange().data()[0];
    auto tr1 = matrix.trange().data()[1];

    auto ncols = tr1.extent();

    if (col_idx < 0 || col_idx >= ncols)
        throw ProgrammingError("Requested col index is out of range", __FILE__, __LINE__);

    auto mat_eig = array_ops::array_to_eigen(matrix);
    Vector<element_type> vec_eig = mat_eig.col(col_idx);

    auto& world = matrix.world();
    TA::TiledRange trange{tr0};
    TA::TensorF norms(trange.tiles_range(), trange.elements_range().volume() / trange.tiles_range().volume());
    TA::SparseShape<float> shape(world, norms, trange);

    TArray result(world, trange, shape);

    auto vec_to_tile = [](TA::Range &range, const Vector<element_type>& arg_vec){
        auto tile = Tile(range);
        const auto p0 = tile.range().lobound()[0];
        const auto pn = tile.range().upbound()[0];

        auto tile_ord = 0;
        for (auto p = p0; p < pn; ++p, ++tile_ord) {
            const auto element = arg_vec[p];
            tile[tile_ord] = element;
        }

        return tile;
    };

    auto const &pmap = result.pmap();
    for (auto it = pmap->begin(); it != pmap->end(); ++it) {
        if (!result.is_zero(*it)) {
            auto range = trange.make_tile_range(*it);
            auto tile = world.taskq.add(vec_to_tile, range, vec_eig);
            result.set(*it, tile);
        }
    }
    world.gop.fence();
    result.truncate();
    return result;
}

}  // namespace detail

}  // namespace gaussian
}  // namespace lcao
}  // namespace mpqc
#endif // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_DECOMPOSED_RIJ_H_
