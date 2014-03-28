//
// cldfgfactory.cpp
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

#include <chemistry/qc/scf/cldfgfactory.hpp>
#include <elemental.hpp>

using namespace mpqc;
using namespace mpqc::TA;
using namespace sc;

using TAMatrix = ClDfGFactory::TAMatrix;
using EMatrix = elem::DistMatrix<double>;

mpqc::TA::ClDfGFactory::ClDfGFactory(sc::Ref<sc::Integral> integral,
                                     sc::Ref<TiledBasisSet> basis,
                                     sc::Ref<TiledBasisSet> dfbasis,
                                     const TAMatrix& density,
                                     sc::Ref<World> world) :
        eri3_(integral->electron_repulsion3()), eri2_(
                integral->electron_repulsion2()), basis_(basis), dfbasis_(
                dfbasis), density_(density), world_(world) {
}

TiledArray::expressions::TensorExpression<TAMatrix::eval_type>
mpqc::TA::ClDfGFactory::operator ()( const std::string& v) {

  if (!df_ints_.is_initialized()) {
    compute_symetric_df_ints();
  }

  auto expr = 2 * (df_ints_("i,j,Q") * (density_("n,m") * df_ints_("n,m,Q")))
                - (df_ints_("i,n,Q") * (density_("n,m") * df_ints_("j,m,Q")));

  return expr;
}

void mpqc::TA::ClDfGFactory::compute_symetric_df_ints() {
  // Get three center ints
  std::shared_ptr<decltype(eri3_)> eri3_ptr(&eri3_);
  // Using the df_ints as temporary storage for the twobody three center ints
  df_ints_ =  Integrals(*world_->madworld(), eri3_ptr, basis_, dfbasis_);

  // Get two center ints
  std::shared_ptr<decltype(eri2_)> eri2_ptr(&eri2_);
  TAMatrix eri2_ints = Integrals(*world_->madworld(), eri2_ptr, dfbasis_);

  // Copy two body two center ints to elemental
  EMatrix eri2_elem = ::TiledArray::array_to_elem(eri2_ints, *world_->elemGrid());

  // Perform the cholesky inverse of the matrix
  elem::HPSDCholesky(elem::LOWER, eri2_elem);
  elem::MakeTriangular(elem::LOWER, eri2_elem);
  elem::TriangularInverse(elem::LOWER, elem::NON_UNIT, eri2_elem);

  // Copy back to TA
  ::TiledArray::elem_to_array(eri2_ints, eri2_elem);

  // Create df_ints this should create a matrix such that X(i,j,Q) = X(j,i,Q)
  // Must check this throughly though
  df_ints_("i,j,Q") = df_ints_("i,j,P") * eri2_ints("P,Q");
}

