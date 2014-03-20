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

#ifndef MPQC_SRC_LIB_UTIL_ELEMENTAL_TA_INTERFACE_HPP
#define MPQC_SRC_LIB_UTIL_ELEMENTAL_TA_INTERFACE_HPP

#include <tiled_array.h>
#include <elemental.hpp>
#include <util/misc/assert.h>

namespace mpqc{

  template <typename T>
  elem::DistMatrix<T> array_to_distmat(TiledArray::Array<T,2> &x,
                                         elem::Grid &g){
    // Get the last tile for size information
    std::vector<unsigned int> sizes = x.trange().elements().size();
    MPQC_ASSERT(sizes.size() == 2);

    // construct DistMat
    elem::DistMatrix<T> mat(sizes[0],sizes[1],g);
    elem::Zero(mat);

    // Create the Axpy interface used to fill DistMat
    elem::AxpyInterface<T> interface;
    // Attach matrix to the interface
    interface.Attach(elem::LOCAL_TO_GLOBAL, mat);

    // Get TA iterator
    auto it = x.begin();
    auto end = x.end();
    for(;it != end; ++it){
      // Get tile and offset info
      typename ::TiledArray::Array<T,2>::value_type tile = *it;
      int t0start = tile.range().start()[0];
      int t1start = tile.range().start()[1];
      int t0size = tile.range().size()[0];
      int t1size = tile.range().size()[1];

      // Create a local elem::Matrix
      elem::Matrix<T> ElemBlock;
      // Attach the tile to it.
      ElemBlock.Attach(t1size, t0size, tile.data(),t0size);

      // Add it to our Distmat
      interface.Axpy(1.0, ElemBlock, t0start, t1start);
    }
    x.get_world().gop.fence();
    interface.Detach();

    return mat;
  }

  template <typename T>
  void distmat_to_array(TiledArray::Array<T,2> &x, elem::DistMatrix<T> &mat){

    elem::AxpyInterface<T> interface;
    interface.Attach(elem::GLOBAL_TO_LOCAL, mat);

    auto it = x.begin();
    auto end = x.end();

    for(;it != end; ++it){
      typename TiledArray::Array<T,2>::value_type tile = *it;
      int t0start = tile.range().start()[0];
      int t1start = tile.range().start()[1];
      int t0size = tile.range().size()[0];
      int t1size = tile.range().size()[1];

      elem::Matrix<double> mat;
      mat.Attach(t1size,t0size,tile.data(),t0size);
      std::fill(mat.Buffer(), mat.Buffer()+t0size*t1size, T(0));

      interface.Axpy(T(1.0), mat, t0start, t1start);
    }
    interface.Detach();
  }



} // namespace mpqc

#endif /* MPQC_SRC_LIB_UTIL_ELEMENTAL_TA_INTERFACE_HPP */
