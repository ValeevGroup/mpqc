//
// orbitalspace_utils.timpl.h
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

#ifdef __GNUG__
#pragma interface
#endif

#ifndef _mpqc_src_lib_chemistry_qc_mbptr12_orbitalspace_utils_timpl_h
#define _mpqc_src_lib_chemistry_qc_mbptr12_orbitalspace_utils_timpl_h

#include <chemistry/qc/basis/obintfactory.h>

namespace sc {

  /// compute one-body integral matrix between in the basis of space
  template <IntegralOneBodyMethod IntMethod>
  RefSymmSCMatrix
  compute_obints(const Ref<OrbitalSpace>& space) {
    Ref<Integral> ints = space->integral()->clone();
    ints->set_basis(space->basis());
    Ref<PetiteList> plist = ints->petite_list();
    RefSymmSCMatrix O_so = compute_onebody_matrix<IntMethod>(plist);
    RefSymmSCMatrix O_ao = plist->to_AO_basis(O_so);
    RefSCMatrix C = space->coefs();
    RefSymmSCMatrix O = O_so.kit()->symmmatrix(C.coldim());
    O.assign(0.0);
    O.accumulate_transform(C, O_ao, SCMatrix::TransposeTransform);
    return O;
  }

  /// compute one-body integral matrix in the space of space_bra and space_ket (space_bra != space_ket)
  template <IntegralOneBodyMethod IntMethod>
  RefSCMatrix
  compute_obints(const Ref<OrbitalSpace>& space_bra, const Ref<OrbitalSpace>& space_ket) {
    throw ProgrammingError("compute_obints not yet implemented",
                           __FILE__, __LINE__);
  }

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
