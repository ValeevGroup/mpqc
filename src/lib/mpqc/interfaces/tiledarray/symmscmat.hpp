//
// symmscmat.hpp
//
// Copyright (C) 2013 Drew Lewis
//
// Authors: Drew Lewis
// Maintainer: Drew Lewis and Edward Valeev
//
// This file is part of the MPQC Toolkit.
//
// The MPQC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The MPQC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the MPQC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifndef mpqc_interfaces_tiledarray_symmscmat_hpp
#define mpqc_interfaces_tiledarray_symmscmat_hpp

#include <tiledarray.h>

namespace mpqc {
/// addtogroup TiledArrayInterface
/// @{

#ifndef DOXYGEN
    namespace detail{
        using TRange1 = ::TiledArray::TiledRange1;
        template<unsigned int N>
        ::TiledArray::Array<double, N>
        init_ta_array(madness::World &world, const TRange1 &trange1){

            std::array<TRange1, N> blocking;
            for(std::size_t i = 0; i < N; ++i){
                blocking[i] = trange1;
            }

            ::TiledArray::TiledRange trange(blocking.begin(), blocking.end());

            ::TiledArray::Array<double, N> array(world, trange);
            return array;
        }

        template<typename RefMat>
        ::TiledArray::Array<double, 2>::value_type
        mat_to_array_task(const ::TiledArray::Range &range, const RefMat &matrix){

            /// This is the tile which we are going to return.
            ::TiledArray::Array<double, 2>::value_type tile(range);

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
        mat_to_array(const RefMat &matrix, ::TiledArray::Array<double, N> &array){

            auto it = array.get_pmap()->begin();
            auto end = array.get_pmap()->end();

            madness::World &world = array.get_world();

            static_assert(N==2, "Feature Not Yet Implimented: Can only copy"
                          " arrays with rank 2");

            for(; it != end; ++it){
                madness::Future<typename ::TiledArray::Array<double,N>::value_type > tile =
                                world.taskq.add(
                                    &mat_to_array_task<RefMat>,
                                    array.trange().make_tile_range(*it),
                                    matrix);

                array.set(*it, tile);

            }

        }
    }// namespace details
#endif // DOXYGEN

    /**
     * Function to copy a sc::SymmSCMatrix to TiledArray::Array
     * @Warning This function is only for testing not for general use
     */
    inline ::TiledArray::Array<double, 2>
    SymmScMat_To_TiledArray(madness::World &world,
                            const sc::Ref<sc::SymmSCMatrix> &matrix,
                            const ::TiledArray::TiledRange1 trange1){

        ::TiledArray::Array<double, 2> ta_array = detail::init_ta_array<2>(world, trange1);

        detail::mat_to_array(matrix, ta_array);

        return ta_array;
    }

    /// @}
    // End TiledArray Interface

} // namespace mpqc









#endif /* mpqc_interfaces_tiledarray_symmscmat_hpp */
