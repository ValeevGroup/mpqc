/*
 * scmat.hpp
 *
 *  Created on: Aug 12, 2013
 *      Author: drewlewis
 */

#ifndef mpqc_interfaces_tiledarray_symmscmat_hpp
#define mpqc_interfaces_tiledarray_symmscmat_hpp

#include <mpqc/interfaces/tiledarray/tiling/trange1.hpp>
#include <tiled_array.h>

namespace TA = TiledArray;
namespace mpqc {

    template<unsigned int N>
    TA::Array<double, N>
    init_ta_array(madness::World &world, const tiling::TRange1 &trange1){

        std::array<tiling::TRange1, N> blocking;
        for(std::size_t i = 0; i < N; ++i){
            blocking[i] = trange1;
        }

        TA::TiledRange trange(blocking.begin(), blocking.end());

        TA::Array<double, N> array(world, trange);
        return array;
    }

    template<typename RefMat>
    TA::Array<double, 2>::value_type
    mat_to_array_task(const TA::Range &range, const RefMat &matrix){

        /// This is the tile which we are going to return.
        TA::Array<double, 2>::value_type tile(range);

        int t0start = tile.range().start()[0];
        int t1start = tile.range().start()[1];
        int t0size = tile.range().size()[0];
        int t1size = tile.range().size()[1];
        int t0end = t0start + t0size;
        int t1end = t1start + t1size;

        for(int i = t0start; i < t0end; ++i){

            int tile_i = i - t0start;
            for(int j = t1start; j < t1end; ++j){
                int tile_j = j - t1start;
                tile[tile_i * t1size + tile_j] = matrix->get_element(i,j);
            }
        }

        return tile;
    }



    template<unsigned int N, typename RefMat>
    void
    mat_to_array(const RefMat &matrix, TA::Array<double, N> &array){

        auto it = array.get_pmap()->begin();
        auto end = array.get_pmap()->end();

        madness::World &world = array.get_world();

        static_assert(N==2, "Feature Not Yet Implimented: Can only copy"
                      " arrays with rank 2");

        for(; it != end; ++it){
            madness::Future<typename TA::Array<double,N>::value_type > tile =
                            world.taskq.add(
                                &mat_to_array_task<RefMat>,
                                array.trange().make_tile_range(*it),
                                matrix);

            array.set(*it, tile);

        }

    }

    TA::Array<double, 2>
    SymmScMat_To_TiledArray(madness::World &world,
                            const sc::Ref<sc::SymmSCMatrix> &matrix,
                            const tiling::TRange1 &trange1){

        TA::Array<double, 2> ta_array = init_ta_array<2>(world, trange1);

        mat_to_array(matrix, ta_array);

        return ta_array;
    }

} // namespace mpqc









#endif /* mpqc_interfaces_tiledarray_symmscmat_hpp */
