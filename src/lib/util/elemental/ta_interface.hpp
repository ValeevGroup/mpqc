//
// ta_interface.hpp
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

#ifdef HAVE_ELEMENTAL
#ifndef MPQC_SRC_LIB_UTIL_ELEMENTAL_TA_INTERFACE_HPP
#define MPQC_SRC_LIB_UTIL_ELEMENTAL_TA_INTERFACE_HPP

#include <tiled_array.h>
#include <elemental.hpp>

namespace mpqc{

  template <typename T>
  elem::DistMatrix<T> array_to_distmat(TiledArray::Array<T,2> x,
                                         elem::Grid &g){
    // Get the last tile for size information
    auto end = x.get_pmap()->end();
    auto tile = x.find(*(end - 1)).get();
    // compute the number of rows and cols
    int n_row = tile.range().start()[0] + tile.range().size()[0];
    int n_col = tile.range().start()[1] + tile.range().size()[1];

    // construct DistMat
    elem::DistMatrix<T> mat(n_row,n_col,g);
    elem::Zero(mat);

    // Create the Axpy interface used to fill DistMat
    elem::AxpyInterface<T> interface;
    // Attach matrix to the interface
    interface.Attach(elem::LOCAL_TO_GLOBAL, mat);

    // Get TA iterator
    auto it = x.get_pmap()->begin();
    for(;it != end; ++it){
      // Get tile and offset info
      auto tile = x.find(*it).get();
      int t0start = tile.range().start()[0];
      int t1start = tile.range().start()[1];
      int t0size = tile.range().size()[0];
      int t1size = tile.range().size()[1];

      // Create a local elem::Matrix
      elem::Matrix<T> ElemBlock;
      // Attach the tile to it.
      ElemBlock.Attach(t0size, t1size, tile.begin(),0);
      // Add it to our Distmat
      interface.Axpy(1.0, ElemBlock, t0start, t1start);
    }
    x.get_world().gop.fence();
    interface.Detach();

    return mat;
  }

} // namespace mpqc

#endif /* MPQC_SRC_LIB_UTIL_ELEMENTAL_TA_INTERFACE_HPP */
#endif // Have_ELEMENTAL
