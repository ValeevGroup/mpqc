//
// r12_amps.h
//
// Copyright (C) 2004 Edward Valeev
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

#ifndef _chemistry_qc_mbptr12_r12amps_h
#define _chemistry_qc_mbptr12_r12amps_h

#ifdef __GNUC__
#pragma interface
#endif

#include <stdexcept>
#include <math/scmat/matrix.h>
#include <chemistry/qc/mbptr12/spin.h>

namespace sc {
  
  class R12IntEval;

/** R12Amplitudes gives the amplitudes of some R12-ansatz-related terms in wave function. The first-order wave function
    terms which result from F12 terms are:
    F<sub>ij</sub><sup>(1)</sup> = C<sub>kl</sub><sup>ij</sup> ( f<sub>12</sub> |kl> - 0.5 f<sub>ab</sub><sup>kl</sup> |ab> - 0.5 f<sub>mn</sub><sup>kl</sup> |mn> - f<sub>am</sub><sup>kl</sup> |am> - r<sub>a'm</sub><sup>kl</sup> |a'm> )
    where f12 is the correlation factor, C are optimal first-order coefficients, and
    f are antisymmetrized integrals over f12 operator. Indices a, b are virtual MOs; m,n are occupied MOs;
    i, j, k, l are active occupied MOs, a' is an RI basis index. */

class R12Amplitudes : public RefCount {
  public:
  R12Amplitudes(const Ref<R12IntEval>& r12eval);
  ~R12Amplitudes();
  
  const RefSCMatrix& T2(SpinCase2 S);
  const RefSCMatrix& Fvv(SpinCase2 S);
  const RefSCMatrix& Foo(SpinCase2 S);
  const RefSCMatrix& Fov(SpinCase2 S);
  const RefSCMatrix& Fox(SpinCase2 S);
  const RefSCMatrix& Fvo(SpinCase2 S);
  const RefSCMatrix& Fxo(SpinCase2 S);

  private:
  Ref<R12IntEval> r12eval_;
  RefSCMatrix T2_[NSpinCases2];
  RefSCMatrix Fvv_[NSpinCases2];
  RefSCMatrix Foo_[NSpinCases2];
  RefSCMatrix Fov_[NSpinCases2];
  RefSCMatrix Fox_[NSpinCases2];
  RefSCMatrix Fvo_[NSpinCases2];
  RefSCMatrix Fxo_[NSpinCases2];
  
  void compute_(SpinCase2 sc2);

};

}
  
#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
