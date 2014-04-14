//
// eigensolver.hpp
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

#ifndef MPQC_MATH_ELEMENTAL_EIGENSOLVER_HPP
#define MPQC_MATH_ELEMENTAL_EIGENSOLVER_HPP

#include <tiledarray_fwd.h>

namespace mpqc {
namespace TA {

  using TiledArray::TArray2D;

  /***
   * Solve generalized Hermitian eigenproblem to form density.
   * This function will only compute N = occupied orbitals number of eigenvalues
   * and vectors.   Thus it can't be used for full MO integral transformation.
   */
  void
  eigensolver_D(const TArray2D &F, const TArray2D &S, TArray2D &D, int occ);

  /***
   * Solve generalized Hermitian eigenproblem for coefficents
   */
  void
  eigensolver_Coeff(const TArray2D &F, const TArray2D &S, TArray2D &C);


} // namespace mpqc::TA
} // namespace mpqc


#endif /* MPQC_MATH_ELEMENTAL_EIGENSOLVER_HPP */
