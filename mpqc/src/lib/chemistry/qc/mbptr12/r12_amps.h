//
// r12_amps.h
//
// Copyright (C) 2004 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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
#include <chemistry/qc/mbptr12/moindexspace.h>

namespace sc {

/** R12Amplitudes gives the amplitudes of some linear-R12-ansatz-related terms in wave function. The first-order wave function
    terms which result from linear R12 terms are:
    F<sub>ij</sub><sup>(1)</sup> = C<sub>kl</sub><sup>ij</sup> ( r<sub>12</sub> |kl> - 0.5 r<sub>ab</sub><sup>kl</sup> |ab> - 0.5 r<sub>mn</sub><sup>kl</sup> |mn> - r<sub>am</sub><sup>kl</sup> |am> - r<sub>a'm</sub><sup>kl</sup> |a'm> )
    where C are optimal first-order coefficients and r are antisymmetrized integrals over r12 operator. Indices a, b are virtual MOs; m,n are occupied MOs;
    i, j, k, l are active occupied MOs, a' is an RI basis index. */

class R12Amplitudes : public RefCount {

  RefSCMatrix T2_aa_, T2_ab_;
  RefSCMatrix Rvv_aa_, Rvv_ab_;
  RefSCMatrix Roo_aa_, Roo_ab_;
  RefSCMatrix Rvo_aa_, Rvo_ab_;
  RefSCMatrix Rxo_aa_, Rxo_ab_;

  public:
    R12Amplitudes(const RefSCMatrix& T2_aa, const RefSCMatrix& T2_ab,
                  const RefSCMatrix& Rvv_aa, const RefSCMatrix& Rvv_ab,
                  const RefSCMatrix& Roo_aa, const RefSCMatrix& Roo_ab,
                  const RefSCMatrix& Rvo_aa, const RefSCMatrix& Rvo_ab,
                  const RefSCMatrix& Rxo_aa, const RefSCMatrix& Rxo_ab)
  {
    T2_aa_ = Rvv_aa;    T2_ab_ = Rvv_ab;
    Rvv_aa_ = Rvv_aa;   Rvv_ab_ = Rvv_ab;
    Roo_aa_ = Roo_aa;   Roo_ab_ = Roo_ab;
    Rvo_aa_ = Rvo_aa;   Rvo_ab_ = Rvo_ab;
    Rxo_aa_ = Rxo_aa;   Rxo_ab_ = Rxo_ab;
  }
    ~R12Amplitudes() {};

    const RefSCMatrix T2_aa() const;
    const RefSCMatrix T2_ab() const;
    const RefSCMatrix Rvv_aa() const;
    const RefSCMatrix Rvv_ab() const;
    const RefSCMatrix Roo_aa() const;
    const RefSCMatrix Roo_ab() const;
    const RefSCMatrix Rvo_aa() const;
    const RefSCMatrix Rvo_ab() const;
    const RefSCMatrix Rxo_aa() const;
    const RefSCMatrix Rxo_ab() const;

};

}
  
#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
