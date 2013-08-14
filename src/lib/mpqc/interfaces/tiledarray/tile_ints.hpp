//
// tile_ints.hpp
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

#ifndef mpqc_interfaces_tiledarray_tileints_h
#define mpqc_interfaces_tiledarray_tileints_h

#include <tiled_array.h>
#include <mpqc/interfaces/tiledarray/trange1.hpp>
#include <chemistry/qc/basis/integral.h>
#include <boost/ref.hpp>
#include <mpqc/utility/foreach.hpp>
#include <mpqc/integrals/integrals.hpp>

namespace mpqc{

    /// Helper class to return number of dimensions of an integral type
    template <typename IntegralEngine>
    struct EngineTypeTraits;
    template <> struct EngineTypeTraits<sc::Ref<sc::TwoBodyInt> > {
        static const size_t ncenters = 4u;
    };
    template <> struct EngineTypeTraits<sc::Ref<sc::TwoBodyThreeCenterInt> > {
        static const size_t ncenters = 3u;
    };
    template <> struct EngineTypeTraits<sc::Ref<sc::TwoBodyTwoCenterInt> > {
        static const size_t ncenters = 2u;
    };
    template <> struct EngineTypeTraits<sc::Ref<sc::OneBodyInt> > {
        static const size_t ncenters = 2u;
    };
    template <> struct EngineTypeTraits<sc::Ref<sc::OneBodyOneCenterInt> > {
        static const size_t ncenters = 1u;
    };

    namespace int_details {
        typedef std::vector<int> shell_range;

        template<typename RefEngine>
        void
        get_integrals(const std::vector<shell_range> &s,
                      RefEngine &engine,
                      TensorRef<double, 2, TensorRowMajor > &tile_map){
            mpqc::integrals::evaluate(engine, s[0], s[1], tile_map);
        }

        template<typename RefEngine>
        void
        get_integrals(const std::vector<shell_range> &s,
                      RefEngine &engine,
                      TensorRef<double, 3, TensorRowMajor > &tile_map){
            mpqc::integrals::evaluate(engine, s[0], s[1], s[2], tile_map);
        }

        template<typename RefEngine>
        void
        get_integrals(const std::vector<shell_range> &s,
                      RefEngine &engine,
                      TensorRef<double, 4, TensorRowMajor > &tile_map){
            mpqc::integrals::evaluate(engine, s[0], s[1], s[2], s[3], tile_map);
        }


    } // namespace int_details

    /**
     * get_integrals takes a tile and a integral engine and calculates the
     * shell offsets for that tile.  It then offloads the filling of the tile
     * to the mpqc::TensorRef and mpqc::Integrals tools.
     */
    template<typename Tile, typename RefEngine>
    void
    get_integrals(Tile &tile, RefEngine &engine){

        // Calculate the rank of our tiles
        constexpr std::size_t rank = EngineTypeTraits<RefEngine>::ncenters;
        typedef int_details::shell_range shell_range;

        std::vector<shell_range> ranges(rank);

        // Loop over each dimension of the tile.
        for(auto i = 0; i < rank; ++i){

            // Get the global indices of the first and last function in a tile
            // This corresponds to basis functions.
            const std::size_t first_func = tile.range().start()[i];
            const std::size_t last_func = tile.range().finish()[i] - 1;

            // Compute the first and  last shell for a given dimension of a tile
            const std::size_t first_shell =
                           engine->basis(i)->function_to_shell(first_func);
            const std::size_t last_shell =
                           engine->basis(i)->function_to_shell(last_func);

            // fill a vector with the shell indices that belong to dimension i
            for(auto j = first_shell; j <  last_shell + 1; ++j){
                ranges[i].push_back(j);
            }
        }

        //  Store the size and dimension of the tile for mpqc::TensorRef
        std::size_t dim[rank];
        //std::copy(tile.range().size().begin(), tile.range().size().end(),
        //          dim);
        for(auto i = 0; i < rank; ++i){
            dim[i] = tile.range().size()[i];
        }

        //const std::size_t (&dim)[rank] = *reinterpret_cast<const std::size_t(*)[rank]>(
        //                        tile.range().size().data());

        TensorRef<double, rank, TensorRowMajor > tile_map(tile.data(), dim);


        int_details::get_integrals(ranges, engine, tile_map);

    }

    /**
     * make_integral_task spawns tasks to fill tiles with integrals.
     */
    template<typename RefPool, typename It, class A>
    void make_integral_task(It first, It last, const A &array,
                       RefPool pool);

    /**
     * integral_task takes a boost reference to a  IntegralEnginePool and
     * first spawns new tasks until the problem is the correct size.
     * Then it gets a thread local instance of an integral engine and uses it
     * to fill tiles with data
     */
    template<typename RefPool, typename It, class A>
    void
    integral_task(It first, It last, A &array,
                  RefPool pool){

        // Divide work.   If the distance between the first and last pointer
        // is greater than some number then split the range in half and send
        // half of the work to another task.  Continue to do this until
        // the distance between the pointers is small enough.
        while(last - first > 5){
            It middle = first;
            std::advance(middle, std::distance(first,last)/2);
            make_integral_task(middle, last, array, pool);
            last = middle;
        }


        // Once the pointer range is small enough unwrap the engine type
        // and get a local instance
        typedef typename boost::unwrap_reference<RefPool>::type PoolType;
        typename PoolType::engine_type engine = pool.get().instance();

        // Loop over the iterator range and  create tiles to populate the
        // TiledArray. Fill the tiles with data in get_integrals
        for(; first != last; ++first){
            typename A::value_type tile(
                            array.trange().make_tile_range(*first));
            get_integrals(tile, engine);

            array.set(*first, tile);
        }
    }

    template<typename RefPool, typename It, class A>
    void
    make_integral_task(It first, It last, const A &array,
                       RefPool pool){
        array.get_world().taskq.add(&integral_task<RefPool, It, A>, first,
                                    last, array, pool);
    }

    /**
     * Initial function called to fill a TiledArray with integrals.
     * It gets the pointers to the first and last tile and then sends the work
     * off to tasks.
     */
    template<typename IntEngPool, class A>
    void fill_tiles(A &array, const IntEngPool &pool){

        // get pointers
        auto begin = array.get_pmap()->begin();
        auto end = array.get_pmap()->end();

        make_integral_task(begin, end, array, boost::cref(pool));
    }

} // namespace mpqc



#endif /* mpqc_interfaces_tiledarray_tileints_h */
