//
// obint.h
//
// Copyright (C) 2001 Edward Valeev
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

#ifndef _chemistry_qc_cints_obint_h
#define _chemistry_qc_cints_obint_h

#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/cints/int1e.h>

namespace sc {

/** This implements most one body integrals in the Cints library. It is
    given a function pointer to the Int1e member that computes the
    particular integral of interest. */
class OneBodyIntCints : public OneBodyInt {
    Ref<Int1eCints> int1ecints_;
    typedef void (Int1eCints::*IntegralFunction)(int,int);
    IntegralFunction intfunc_;
  public:
    OneBodyIntCints(Integral*,
                 const Ref<GaussianBasisSet>&, const Ref<GaussianBasisSet>&,
                 IntegralFunction);
    ~OneBodyIntCints();
 
    void set_multipole_origin(const Ref<DipoleData>&);
    void set_EdotV_origin(const Ref<EfieldDotVectorData>&);
    void set_Q_origin(const Ref<PointChargeData>&);

    void compute_shell(int,int);

    bool cloneable() const;
    Ref<OneBodyInt> clone();
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
