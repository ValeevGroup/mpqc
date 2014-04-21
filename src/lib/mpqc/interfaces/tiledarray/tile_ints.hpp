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

#ifndef mpqc_interfaces_tiledarray_tileints_hpp
#define mpqc_interfaces_tiledarray_tileints_hpp

#include <tiledarray.h>
#include <mpqc/interfaces/tiledarray/tiling/trange1.hpp>
#include <chemistry/qc/basis/integral.h>
#include <boost/ref.hpp>
#include <mpqc/utility/foreach.hpp>
#include <mpqc/integrals/integrals.hpp>

namespace mpqc{

    /// @addtogroup ChemistryBasisIntegralTA
    /// @{

#ifndef DOXYGEN
    // Helper class to get Integral Engine typetraits.
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

    /*
     * Does the dirty work of getting integrals into TA::Tiles
     */
    namespace int_details {
        typedef std::vector<int> shell_range;

        // Takes a shell_range, sc::Engine_Type, and a TensorRef map to a TA::Tile
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


    /*
     * Computes shell ranges to pass to a TensorRef for getting integrals
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

        // passes the TiledArray size into a fixed c-style array for Tensor class
        const std::size_t (&dim)[rank] =
                        *reinterpret_cast<const std::size_t(*)[rank]>(
                                            tile.range().size().data());

        TensorRef<double, rank, TensorRowMajor > tile_map(tile.data(), dim);


        int_details::get_integrals(ranges, engine, tile_map);

    }

    template<typename RefPool, typename It, class A>
    void make_integral_task(It first, It last, const A &array,
                       RefPool pool);

    /*
     * Splits work in into manageable chunks by creating tasks. And then fills
     * each tile with integrals.
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

    /*
     * Spawns tasks to fill tiles with integrals.
     */
    template<typename RefPool, typename It, class A>
    void
    make_integral_task(It first, It last, const A &array,
                       RefPool pool){
        array.get_world().taskq.add(&integral_task<RefPool, It, A>, first,
                                    last, array, pool);
    }

#endif //DOXYGEN

    /**
     * Initial function called to fill a TiledArray with integrals.
     * @param[in,out] array is a TiledArray::Array that will be filled with data
     * @param[in] pool is an IntegralEnginePool object to provide integrals.
     */
    template<typename IntEngPool, class A>
    void fill_tiles(A &array, const IntEngPool &pool){

        // get pointers
        auto begin = array.get_pmap()->begin();
        auto end = array.get_pmap()->end();

        // Create tasks to fill tiles with data. Boost const reference is used
        // because Integral Engine pool is not copyable, but when sent to the
        // Madness task queue all objects are copied.
        make_integral_task(begin, end, array, boost::cref(pool));
    }

} // namespace mpqc


/// @}

#endif /* mpqc_interfaces_tiledarray_tileints_hpp */
