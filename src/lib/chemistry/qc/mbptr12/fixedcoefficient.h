//
// fixedcoefficient.h
//
// Copyright (C) 2008 Martin Torheyden
//
// Author: Martin Torheyden <mtorhey@vt.edu>
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

#ifndef _chemistry_qc_mbptr12_fixedcoefficient_h_
#define _chemistry_qc_mbptr12_fixedcoefficient_h_

#include <util/ref/ref.h>
#include <math/scmat/matrix.h>
#include <util/class/scexception.h>
#include <chemistry/qc/wfn/spin.h>
#include <chemistry/qc/mbptr12/r12technology.h>
#include <math/mmisc/pairiter.h>

namespace sc {

/**
 * Computes fixed coefficients determined
 * according to the cusp conditions for
 * geminal (r12-dependent) functions that have been
 * normalized so that coefficient of r12 in Taylor expansion around r12=0 is 1
 */
class CuspConsistentGeminalCoefficient : virtual public RefCount {
  private:
    SpinCase2 pairspin_;
  public:
    CuspConsistentGeminalCoefficient(SpinCase2 pairspin);
    ~CuspConsistentGeminalCoefficient(){}
    /**
     * OW: Indices of geminal generating space
     * PQ: Indices of space from which geminal substitutions are allowed
     */
    double C(unsigned int O, unsigned int W,
             unsigned int P, unsigned int Q);
};

}

#endif /*_chemistry_qc_mbptr12_fixedcoefficient_h_*/
