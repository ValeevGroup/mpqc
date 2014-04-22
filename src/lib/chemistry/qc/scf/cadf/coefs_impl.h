//
// coefs_impl.h
//
// Copyright (C) 2014 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: Apr 16, 2014
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

#ifndef _chemistry_qc_scf_coefs_impl_h
#define _chemistry_qc_scf_coefs_impl_h

#include "iters.h"
#include "cadfclhf.h"

namespace sc {

template<template<typename...> class container, typename map_type>
void
CADFCLHF::get_coefs_ish_jsh(
    const ShellData& ish,
    const ShellData& jsh,
    int ithr,
    container<map_type>& coefsA,
    container<map_type>& coefsB
) {
  const int dfnbfAB = ish.center == jsh.center ? ish.atom_dfnbf : ish.atom_dfnbf + jsh.atom_dfnbf;

  // Probably more efficient to just do all of the decompositions
  //   in threads beforehand, but oh well
  std::shared_ptr<Decomposition> decomp = get_decomposition(                       //latex `\label{sc:coefgetdecomp}`
      ish, jsh, metric_ints_2c_[ithr]
  );

  Eigen::MatrixXd ij_M_X(ish.nbf*jsh.nbf, dfnbfAB);
  for(auto&& ksh : iter_shells_on_center(dfbs_, ish.center)){
    auto ij_M_k = ints_to_eigen(
        ish, jsh, ksh,
        metric_ints_3c_[ithr],
        metric_oper_type_
    );
    // Surely this transfer can be done more efficiently
    for(auto&& ibf : function_range(ish)){
      for(auto&& jbf : function_range(jsh)){
        const int ijbf = ibf.bfoff_in_shell * jsh.nbf + jbf.bfoff_in_shell;
        ij_M_X.row(ijbf).segment(ksh.bfoff_in_atom, ksh.nbf) = ij_M_k->row(ijbf);
      } // end loop over functions in jsh
    } // end loop over functions in ish
  } // end loop over shells on ish.center

  if(ish.center != jsh.center){
    for(auto&& ksh : iter_shells_on_center(dfbs_, jsh.center)){
      auto ij_M_k = ints_to_eigen(
          ish, jsh, ksh,
          metric_ints_3c_[ithr],
          metric_oper_type_
      );
      for(auto&& ibf : function_range(ish)){
        for(auto&& jbf : function_range(jsh)){
          const int ijbf = ibf.bfoff_in_shell * jsh.nbf + jbf.bfoff_in_shell;
          const int dfbfoff = ish.atom_dfnbf + ksh.bfoff_in_atom;
          ij_M_X.row(ijbf).segment(dfbfoff, ksh.nbf) = ij_M_k->row(ijbf);
        } // end loop over functions in jsh
      } // end loop over functions in ish
    } // end loop over shells on jsh.center
  } // end if ish.center != jsh.center

  //----------------------------------------//

  //auto& Ctmp = decomp->solve(ij_M_X.transpose());
  //int ijbf = 0;
  //for(auto& Ca : coefsA) {
  //  Ca = Ctmp.col(ijbf++).head(ish.atom_dfnbf);
  //}
  //ijbf = 0;
  //for(auto& Cb : coefsB) {
  //  Cb = Ctmp.col(ijbf++).tail(jsh.atom_dfnbf);
  //}

  for(int ijbf = 0; ijbf < ij_M_X.rows(); ++ijbf) {
    auto& Ctmp = decomp->solve(ij_M_X.row(ijbf).transpose());
    coefsA[ijbf] = Ctmp.head(ish.atom_dfnbf);
    if(!coefsB.empty()) {
      coefsB[ijbf] = Ctmp.tail(jsh.atom_dfnbf);
    }
  }

}

} // end namespace sc

#endif /* _chemistry_qc_scf_coefs_impl_h */
