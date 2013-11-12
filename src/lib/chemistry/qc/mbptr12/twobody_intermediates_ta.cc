//
// twobody_intermediates_ta.cc
//
// Copyright (C) 2013 Edward Valeev
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
#pragma implementation
#endif

#include <cassert>
#include <mpqc_config.h>

#if defined(HAVE_MPQC3_RUNTIME)

#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/sr_r12intermediates.h>

using namespace sc;
namespace TA = TiledArray;

void
R12IntEval::V_diag_ta() {
  SingleReference_R12Intermediates<double> srr12intrmds(madness::World::get_default(),
                                                        this->r12world());
#if 0
  auto Vpair = srr12intrmds.V_diag();
  ExEnv::out0() << indent << "V_ij_ij" << std::endl << Vpair.first << std::endl
                << indent << "V_ij_ji" << std::endl << Vpair.second << std::endl;

  auto Xpair = srr12intrmds.X_diag();
  ExEnv::out0() << indent << "X_ij_ij" << std::endl << Xpair.first << std::endl
                << indent << "X_ij_ji" << std::endl << Xpair.second << std::endl;

  auto Bpair = srr12intrmds.B_diag();
  ExEnv::out0() << indent << "B_ij_ij" << std::endl << Bpair.first << std::endl
                << indent << "B_ij_ji" << std::endl << Bpair.second << std::endl;
#endif

  bool vir_cabs_coupling = true; // need CABS singles into vir+CABS? set to true
  this->compute_emp2_cabs_singles_noncanonical(vir_cabs_coupling);
  srr12intrmds.set_T1_cabs(this->T1_cabs_[Alpha]);

  assert(this->orbital_registry()->key_exists("A'"));

  auto rdm1 = srr12intrmds.rdm1();
}

void
R12IntEval::compute_ccr12_1rdm(const RefSCMatrix& T1, const Ref<DistArray4> (&T2)[NSpinCases2])
{
  SingleReference_R12Intermediates<double> srr12intrmds(madness::World::get_default(),
                                                        this->r12world());
  srr12intrmds.set_T1(T1);
  srr12intrmds.set_T2(T2);
  auto rdm1 = srr12intrmds.rdm1();
}

#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
