//
// orbitalspace_utils.h
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

#ifndef _mpqc_src_lib_chemistry_qc_mbptr12_orbitalspaceutils_h
#define _mpqc_src_lib_chemistry_qc_mbptr12_orbitalspaceutils_h

#include <chemistry/qc/wfn/orbitalspace.h>
#include <chemistry/qc/basis/orthog.h>
#include <chemistry/qc/lcao/fockbuild_runtime.h>
#include <chemistry/qc/basis/obintfactory.h>

namespace sc {

  /** Compute span of bs and create corresponding mospace referred to by name.
      Setting nlindep on input to a non-negative value
      will cause the number of linear dependencies to be exactly its input value,
      if possible (lindep_tol is then ignored).
      This is not implemented for Gram-Schmidt orthogonalization and thus in that case nlindep
      is ignored on input. On output nlindep contains the number of linear dependencies.
      */
   Ref<OrbitalSpace> orthogonalize(const std::string& id, const std::string& name, const Ref<GaussianBasisSet>& bs,
                                   const Ref<Integral>& integral, OverlapOrthog::OrthogMethod orthog_method, double lindep_tol,
                                   int& nlindep);

  /** Project space1 on space2. This routine computes X2 such that C1.S12.X2 = I,
      where I is identity matrix, C1 is space1, and X2 spans
      subspace of space2. X2 is returned. */
  Ref<OrbitalSpace> gen_project(const Ref<OrbitalSpace>& space1, const Ref<OrbitalSpace>& space2,
                                const std::string& id, const std::string& name, double lindep_tol);

  /** Compute subspace X2 of space2 which is orthogonal complement to space1, i.e.,
      C1.S12.X2=0, where 0 is the null matrix.
  */
  Ref<OrbitalSpace> orthog_comp(const Ref<OrbitalSpace>& space1, const Ref<OrbitalSpace>& space2,
                                const std::string& id, const std::string& name, double lindep_tol);

  /// compute one-body integral matrix in the space of space_bra and space_ket (space_bra != space_ket)
  template <IntegralOneBodyMethod IntMethod>
  RefSCMatrix
  compute_obints(const Ref<OrbitalSpace>& space_bra, const Ref<OrbitalSpace>& space_ket);

  /// compute one-body integral matrix between in the basis of space
  template <IntegralOneBodyMethod IntMethod>
  RefSymmSCMatrix
  compute_obints(const Ref<OrbitalSpace>& space);

  /// compute overlap between space1 and space2
  RefSCMatrix
  compute_overlap_ints(const Ref<OrbitalSpace>& space1, const Ref<OrbitalSpace>& space2);

  /// Compute electric dipole and quadrupole moment matrices in the basis of space1 and space2
  void compute_multipole_ints(const Ref<OrbitalSpace>& space1,
                              const Ref<OrbitalSpace>& space2,
                              RefSCMatrix& MX,
                              RefSCMatrix& MY,
                              RefSCMatrix& MZ,
                              RefSCMatrix& MXX,
                              RefSCMatrix& MYY,
                              RefSCMatrix& MZZ,
                              RefSCMatrix& MXY,
                              RefSCMatrix& MXZ,
                              RefSCMatrix& MYZ);

  /// canonicalize A
  Ref<OrbitalSpace>
  compute_canonvir_space(const Ref<FockBuildRuntime>& fb_rtime,
                         const Ref<OrbitalSpace>& A,
                         SpinCase1 spin);

} // end of namespace sc

#include <chemistry/qc/wfn/orbitalspace_utils.timpl.h>

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
