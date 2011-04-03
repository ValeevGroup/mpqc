//
// obintfactory.h
//
// Copyright (C) 2009 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
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

#ifndef _mpqc_src_lib_chemistry_qc_basis_obintfactory_h
#define _mpqc_src_lib_chemistry_qc_basis_obintfactory_h

#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/basis/symmint.h>

namespace sc {

  typedef Ref<OneBodyInt> (Integral::*IntegralOneBodyMethod)();
  /// Creates matrix representation of a one-body operator using point-group symmetry
  /// @tparam IntMethod member function of Integral that produces the evaluator for the desired integrals
  /// @param plist PetiteList object initialized with the proper Integral factory and basis set
  template <IntegralOneBodyMethod IntMethod>
  RefSymmSCMatrix compute_onebody_matrix(Ref<PetiteList>& plist) {
    Ref<GaussianBasisSet> basis = plist->basis();
    Integral* integral = plist->integral().pointer();
    RefSymmSCMatrix O_skel(plist->AO_basisdim(), basis->matrixkit());
    O_skel.assign(0.0);
    Ref<SCElementOp> intop =
        new OneBodyIntOp(new SymmOneBodyIntIter((integral ->* IntMethod)(), plist));
    O_skel.element_op(intop);
    intop = 0;
    RefSymmSCMatrix O(plist->SO_basisdim(), basis->so_matrixkit());
    plist->symmetrize(O_skel, O);
    return O;
  }

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
