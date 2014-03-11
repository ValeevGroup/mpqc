//
// array_ints.hpp
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
#ifndef mpqc_interfaces_tiledarray_integrals_hpp
#define mpqc_interfaces_tiledarray_integrals_hpp

#include "tile_ints.hpp"
#include <chemistry/qc/basis/tiledbasisset.hpp>

namespace mpqc {

  /// @addtogroup ChemistryBasisIntegralTA
  /// @{

  template<typename ShrPtrPool>
  ::TiledArray::Array<double,
      EngineTypeTraits<typename std::pointer_traits<ShrPtrPool>::element_type::engine_type>::ncenters>
  Integrals(
      madness::World &world, const ShrPtrPool &pool,
      const sc::Ref<mpqc::TA::TiledBasisSet> &tbasis) {

    namespace TA = ::TiledArray;

    // Get the the type of integral that we are computing.
    typedef typename std::pointer_traits<ShrPtrPool>::element_type::engine_type engine_type;
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

/// @} // ChemistryBasisIntegralTA

}// namespace mpqc

#endif /* mpqc_interfaces_tiledarray_integrals_hpp */
