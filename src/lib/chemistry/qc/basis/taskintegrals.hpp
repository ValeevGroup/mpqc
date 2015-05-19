//
// taskintegrals.hpp
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

#ifndef MPQC_CHEMISTRY_QC_BASIS_TASKINTEGRALS_HPP_
#define MPQC_CHEMISTRY_QC_BASIS_TASKINTEGRALS_HPP_

#include <tiledarray.h>
#include <memory>
#include <chemistry/qc/basis/tiledbasisset.hpp>
#include <mpqc/integrals/integrals.hpp>
#include <chemistry/qc/basis/integral.h>

namespace mpqc {
  namespace TA {

/// @addtogroup ChemistryBasisIntegralTA
/// @{
    template<typename T>
    using PoolPtrType = typename std::pointer_traits<T>::element_type;

#ifndef DOXYGEN
// Helper class to get Integral Engine typetraits.
    template<typename IntegralEngine>
    struct EngineTypeTraits;
    template<> struct EngineTypeTraits<sc::Ref<sc::TwoBodyInt> > {
      static const size_t ncenters = 4u;
    };
    template<> struct EngineTypeTraits<sc::Ref<sc::TwoBodyThreeCenterInt> > {
      static const size_t ncenters = 3u;
    };
    template<> struct EngineTypeTraits<sc::Ref<sc::TwoBodyTwoCenterInt> > {
      static const size_t ncenters = 2u;
    };
    template<> struct EngineTypeTraits<sc::Ref<sc::OneBodyInt> > {
      static const size_t ncenters = 2u;
    };
    template<> struct EngineTypeTraits<sc::Ref<sc::OneBodyOneCenterInt> > {
      static const size_t ncenters = 1u;
    };

    /*
     * Does the dirty work of getting integrals into TA::Tiles
     */
    namespace int_details {
      typedef std::vector<int> shell_range;

// Takes a shell_range, sc::Engine_Type, and a TensorRef map to a TA::Tile
      template<typename RefEngine>
      void get_integrals(const std::vector<shell_range> &s, RefEngine &engine,
                         TensorRef<double, 2, TensorRowMajor> &tile_map) {
        mpqc::integrals::evaluate(engine, s[0], s[1], tile_map);
      }

      template<typename RefEngine>
      void get_integrals(const std::vector<shell_range> &s, RefEngine &engine,
                         TensorRef<double, 3, TensorRowMajor> &tile_map) {
        mpqc::integrals::evaluate(engine, s[0], s[1], s[2], tile_map);
      }

      template<typename RefEngine>
      void get_integrals(const std::vector<shell_range> &s, RefEngine &engine,
                         TensorRef<double, 4, TensorRowMajor> &tile_map) {
        mpqc::integrals::evaluate(engine, s[0], s[1], s[2], s[3], tile_map);
      }

    } // namespace int_details

    /*
     * Computes shell ranges to pass to a TensorRef for getting integrals
     */
    template<typename Tile, typename RefEngine>
    void get_integrals(Tile &tile, RefEngine &engine) {

      // Calculate the rank of our tiles
      constexpr std::size_t rank = EngineTypeTraits<RefEngine>::ncenters;
      typedef int_details::shell_range shell_range;

      std::vector<shell_range> ranges(rank);

      // Loop over each dimension of the tile.
      for (std::size_t i = 0; i < rank; ++i) {

        // Get the global indices of the first and last function in a tile
        // This corresponds to basis functions.
        const std::size_t first_func = tile.range().start()[i];
        const std::size_t last_func = tile.range().finish()[i] - 1;

        // Compute the first and  last shell for a given dimension of a tile
        const std::size_t first_shell = engine->basis(i)->function_to_shell(
                first_func);
        const std::size_t last_shell = engine->basis(i)->function_to_shell(
                last_func);

        // fill a vector with the shell indices that belong to dimension i
        for (auto j = first_shell; j < last_shell + 1; ++j) {
          ranges[i].push_back(j);
        }
      }

      // passes the TiledArray size into a fixed c-style array for Tensor class
      const std::size_t (&dim)[rank] =
              *reinterpret_cast<const std::size_t (*)[rank]>(tile.range().size());

      TensorRef<double, rank, TensorRowMajor> tile_map(tile.data(), dim);

      int_details::get_integrals(ranges, engine, tile_map);

    }

    template<typename ShrPtrPool, typename It, class A>
    void make_integral_task(It first, It last, const A &array, ShrPtrPool pool);

    /*
     * Splits work in into manageable chunks by creating tasks. And then fills
     * each tile with integrals.
     */
    template<typename ShrPtrPool, typename It, class A>
    void integral_task(It first, It last, A &array, ShrPtrPool &pool) {

      // Divide work.   If the distance between the first and last pointer
      // is greater than some number then split the range in half and send
      // half of the work to another task.  Continue to do this until
      // the distance between the pointers is small enough.
      while (last - first > 10) {
        It middle = first;
        std::advance(middle, std::distance(first, last) / 2);
        make_integral_task(middle, last, array, pool);
        last = middle;
      }

      // Once the pointer range is small enough unwrap the engine type
      // and get a local instance
      typename PoolPtrType<ShrPtrPool>::engine_type engine = pool->instance();

      // Loop over the iterator range and  create tiles to populate the
      // TiledArray. Fill the tiles with data in get_integrals
      for (; first != last; ++first) {
        typename A::value_type tile(array.trange().make_tile_range(*first));
        get_integrals(tile, engine);

        array.set(*first, tile);
      }
    }

    /*
     * Spawns tasks to fill tiles with integrals.
     */
    template<typename ShrPtrPool, typename It, class A>
    void make_integral_task(It first, It last, const A &array, ShrPtrPool pool){
      array.get_world().taskq.add(&integral_task<ShrPtrPool, It, A>, first, last,
                                  array, pool);
    }

#endif //DOXYGEN
    /**
     * Initial function called to fill a TiledArray with integrals.
     * @param[in,out] array is a TiledArray::Array that will be filled with data
     * @param[in] pool is an IntegralEnginePool object to provide integrals.
     */
    template<typename ShrPtrPool, class A>
    void fill_tiles(A &array, const ShrPtrPool &pool) {

      // get pointers
      auto begin = array.get_pmap()->begin();
      auto end = array.get_pmap()->end();

      // Create tasks to fill tiles with data. Boost const reference is used
      // because Integral Engine pool is not copyable, but when sent to the
      // Madness task queue all objects are copied.
      make_integral_task(begin, end, array, pool);
    }

/// @addtogroup ChemistryBasisIntegralTA
/// @{

    template<typename ShrPtrPool>
    ::TiledArray::Array<double,
            EngineTypeTraits<typename PoolPtrType<ShrPtrPool>::engine_type>::ncenters >
    Integrals(
            madness::World &world, const ShrPtrPool &pool,
            const sc::Ref<mpqc::TA::TiledBasisSet> &tbasis) {

      namespace TA = ::TiledArray;

      // Get the the type of integral that we are computing.
      typedef typename PoolPtrType<ShrPtrPool>::engine_type engine_type;
      // Determine the dimensions of our integrals as well as our TiledArray
      constexpr size_t rank = EngineTypeTraits<engine_type>::ncenters;

      // Get the array to initialize the TiledArray::TiledRange using the
      // TiledBasis
      std::array<TiledArray::TiledRange1, rank> blocking;
      for (auto i = 0; i < rank; ++i) {
        blocking[i] = tbasis->trange1();
      }

      // Construct the TiledArray::TiledRange object
      TA::TiledRange trange(blocking.begin(), blocking.end());

      // Initialize the TiledArray
      TA::Array<double, rank> array(world, trange);

      // Fill the TiledArray with data by looping over tiles and sending
      // each tile to a madness task to be filled in parallel.
      fill_tiles(array, pool);

      return array;
    }

    template<typename ShrPtrPool>
    ::TiledArray::Array<double,
            EngineTypeTraits<typename PoolPtrType<ShrPtrPool>::engine_type>::ncenters >
    Integrals( madness::World &world, const ShrPtrPool &pool,
            const sc::Ref<mpqc::TA::TiledBasisSet> &tbasis,
            const sc::Ref<mpqc::TA::TiledBasisSet> &dftbasis) {

      namespace TA = ::TiledArray;

      // Get the the type of integral that we are computing.
      typedef typename PoolPtrType<ShrPtrPool>::engine_type engine_type;
      // Determine the dimensions of our integrals as well as our TiledArray
      constexpr size_t rank = EngineTypeTraits<engine_type>::ncenters;

      // Get the array to initialize the TiledArray::TiledRange using the
      // TiledBasis
      std::array<TiledArray::TiledRange1, rank> blocking;
      for (auto i = 0; i < (rank - 1); ++i) {
        blocking[i] = tbasis->trange1(); // Asign first 2 dims to regular basis
      }
      blocking.back() = dftbasis->trange1(); // Asign last dim to df basis

      // Construct the TiledArray::TiledRange object
      TA::TiledRange trange(blocking.begin(), blocking.end());

      // Initialize the TiledArray
      TA::Array<double, rank> array(world, trange);

      // Fill the TiledArray with data by looping over tiles and sending
      // each tile to a madness task to be filled in parallel.
      fill_tiles(array, pool);

      return array;
    }

/// @} // ChemistryBasisIntegralTA

  }// namespace TA
} // namespace mpqc

#endif /* MPQC_CHEMISTRY_QC_BASIS_TASKINTEGRALS_HPP_ */
