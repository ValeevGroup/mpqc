//
// decomp.cc
//
// Copyright (C) 2014 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: Feb 13, 2014
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//


#include <util/container/conc_cache.h>

#include "cadfclhf.h"

using namespace sc;

shared_ptr<CADFCLHF::Decomposition>
CADFCLHF::get_decomposition(int ish, int jsh, Ref<TwoBodyTwoCenterInt> ints)
{
  // TODO atom swap symmetry?
  const int atomA = basis()->shell_to_center(ish);
  const int atomB = basis()->shell_to_center(jsh);
  const auto& g2 = *g2_full_ptr_;
  //----------------------------------------//
  return decomps_->get(atomA, atomB, [&](){
    // Make the decomposition
    //----------------------------------------//
    const int dfnshA = dfbs_->nshell_on_center(atomA);
    const int dfnbfA = dfbs_->nbasis_on_center(atomA);
    const int dfshoffA = dfbs_->shell_on_center(atomA, 0);
    const int dfbfoffA = dfbs_->shell_to_function(dfshoffA);
    //----------------------------------------//
    const int dfnshB = dfbs_->nshell_on_center(atomB);
    const int dfnbfB = dfbs_->nbasis_on_center(atomB);
    const int dfshoffB = dfbs_->shell_on_center(atomB, 0);
    const int dfbfoffB = dfbs_->shell_to_function(dfshoffB);
    //----------------------------------------//
    // Compute the integrals we need
    const int nrows = atomA == atomB ? dfnbfA : dfnbfA + dfnbfB;
    Eigen::MatrixXd g2AB(nrows, nrows);
    g2AB.block(0, 0, dfnbfA, dfnbfA) = g2.block(
        dfbfoffA, dfbfoffA,
        dfnbfA, dfnbfA
    );
    if(atomA != atomB) {
      g2AB.block(0, dfnbfA, dfnbfA, dfnbfB) = g2.block(
          dfbfoffA, dfbfoffB,
          dfnbfA, dfnbfB
      );
      g2AB.block(dfnbfA, 0, dfnbfB, dfnbfA) = g2.block(
          dfbfoffB, dfbfoffA,
          dfnbfB, dfnbfA
      );
      g2AB.block(dfnbfA, dfnbfA, dfnbfB, dfnbfB) = g2.block(
          dfbfoffB, dfbfoffB,
          dfnbfB, dfnbfB
      );
    }
    return make_shared<Decomposition>(g2AB);
  });
}

//////////////////////////////////////////////////////////////////////////////////

