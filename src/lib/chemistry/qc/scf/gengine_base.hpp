//
// gfactory.hpp
//
// Copyright (C) 2014 Drew Lewis
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

#ifndef MPQC_CHEMISTRY_QC_SCF_GENGINE_HPP
#define MPQC_CHEMISTRY_QC_SCF_GENGINE_HPP

#include <tiledarray_fwd.h>
#include <vector>
#include <util/class/class.h>

// Expression Foward Declaration
namespace TiledArray {
  namespace expressions {
    template <typename>
    class TensorExpression;
  } // namespace expressions
} // namespace TiledArray

namespace mpqc {
  namespace TA {

    class GEngineBase : virtual public sc:: DescribedClass {
    public:
      typedef TiledArray::TArray2D TAMatrix;
      typedef TiledArray::TensorD TATensor;
      typedef TiledArray::expressions::TensorExpression<TATensor> return_type;

      GEngineBase() = default;
      virtual ~GEngineBase() = default;

      virtual // Ensure that the user can set the densities
      void
      set_densities(std::vector<TAMatrix*>) = 0;

      virtual // Return true if the density has been set
      bool
      densities_set() = 0;

      virtual
      return_type
      operator()(const std::string &v) = 0;

    private:
      static sc::ClassDesc class_desc_;

    }; // class GEngine

  } // namespace TA
} // namespace mpqc

#endif /* MPQC_CHEMISTRY_QC_SCF_GENGINE_HPP */
