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
  //----------------------------------------//
  return decomps_->get(atomA, atomB, [&](){
    // Make the decomposition
    std::shared_ptr<Decomposition> decompAB;
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
    // AA integrals
    for(int ishA = dfshoffA; ishA < dfshoffA + dfnshA; ++ishA){
      const int dfbfoffiA = dfbs_->shell_to_function(ishA);
      const int dfnbfiA = dfbs_->shell(ishA).nfunction();
      for(int jshA = dfshoffA; jshA < dfshoffA + dfnshA; ++jshA){
        const int dfbfoffjA = dfbs_->shell_to_function(jshA);
        const int dfnbfjA = dfbs_->shell(jshA).nfunction();
        auto shell_ints = ints_to_eigen(
            ishA, jshA, ints,
            metric_oper_type_
        );
        g2AB.block(
            dfbfoffiA - dfbfoffA, dfbfoffjA - dfbfoffA,
            dfnbfiA, dfnbfjA
        ) = *shell_ints;
        g2AB.block(
            dfbfoffjA - dfbfoffA, dfbfoffiA - dfbfoffA,
            dfnbfjA, dfnbfiA
        ) = shell_ints->transpose();
      }
    }
    if(atomA != atomB) {
      // AB integrals
      for(int ishA = dfshoffA; ishA < dfshoffA + dfnshA; ++ishA){
        const int dfbfoffiA = dfbs_->shell_to_function(ishA);
        const int dfnbfiA = dfbs_->shell(ishA).nfunction();
        for(int jshB = dfshoffB; jshB < dfshoffB + dfnshB; ++jshB){
          const int dfbfoffjB = dfbs_->shell_to_function(jshB);
          const int dfnbfjB = dfbs_->shell(jshB).nfunction();
          auto shell_ints = ints_to_eigen(
              ishA, jshB, ints,
              metric_oper_type_
          );
          g2AB.block(
              dfbfoffiA - dfbfoffA, dfbfoffjB - dfbfoffB + dfnbfA,
              dfnbfiA, dfnbfjB
          ) = *shell_ints;
          g2AB.block(
              dfbfoffjB - dfbfoffB + dfnbfA, dfbfoffiA - dfbfoffA,
              dfnbfjB, dfnbfiA
          ) = shell_ints->transpose();
        }
      }
      // BB integrals
      for(int ishB = dfshoffB; ishB < dfshoffB + dfnshB; ++ishB){
        const int dfbfoffiB = dfbs_->shell_to_function(ishB);
        const int dfnbfiB = dfbs_->shell(ishB).nfunction();
        for(int jshB = dfshoffB; jshB < dfshoffB + dfnshB; ++jshB){
          const int dfbfoffjB = dfbs_->shell_to_function(jshB);
          const int dfnbfjB = dfbs_->shell(jshB).nfunction();
          auto shell_ints = ints_to_eigen(
              ishB, jshB, ints,
              metric_oper_type_
          );
          g2AB.block(
              dfbfoffiB - dfbfoffB + dfnbfA, dfbfoffjB - dfbfoffB + dfnbfA,
              dfnbfiB, dfnbfjB
          ) = *shell_ints;
          g2AB.block(
              dfbfoffjB - dfbfoffB + dfnbfA, dfbfoffiB - dfbfoffB + dfnbfA,
              dfnbfjB, dfnbfiB
          ) = shell_ints->transpose();
        }
      }
    }
    return make_shared<Decomposition>(g2AB);
  });
}

//////////////////////////////////////////////////////////////////////////////////

