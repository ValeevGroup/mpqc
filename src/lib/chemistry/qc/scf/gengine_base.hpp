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

//#include <tiledarray_fwd.h> can't use until fwd gets fixed if possible
#include <tiledarray.h>
#include <vector>
#include <utility>
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
      typedef TiledArray::Array<double,2> TAMatrix;

      GEngineBase() = default;
      virtual ~GEngineBase() = default;

      virtual // Ensure that the user can set the densities
      void
      set_densities(std::vector<TAMatrix*>) = 0;

      virtual // Return true if the densities have been set
      bool
      densities_set() = 0;

      virtual
      void
      set_coefficients(std::vector<TAMatrix*>) = 0;

      virtual // Return true if the coefficients have been set
      bool
      coefficients_set() = 0;

      // If using coefficients for contractions
      virtual
      bool
      using_coeff() = 0;

      // If using densities for contraction
      virtual
      bool
      using_density() = 0;

      virtual
      TAMatrix
      operator()(const std::string v) = 0;

    private:
      static sc::ClassDesc class_desc_;

    }; // class GEngine

  } // namespace TA
} // namespace mpqc

#endif /* MPQC_CHEMISTRY_QC_SCF_GENGINE_HPP */
