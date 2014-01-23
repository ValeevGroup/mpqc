//
// fixedcoefficient.cc
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


#include <chemistry/qc/mbptr12/fixedcoefficient.h>

namespace sc {
  
  CuspConsistentGeminalCoefficient::CuspConsistentGeminalCoefficient(SpinCase2 pairspin)
    : pairspin_(pairspin) { }
  
  double CuspConsistentGeminalCoefficient::C(unsigned int O, unsigned int W,
                                             unsigned int P, unsigned int Q) {
    // only C_ij_ij or C_ij_ji are nonzero
    if( !(((O==P) && (W==Q)) ||
          ((O==Q) && (W==P))) ) {
      return(0.0);
    }

    // cusp conditions are the same depend on spin
    const double Cp_ij_ij = 1.0/2.0;
    const double Cm_ij_ij = 1.0/4.0;
    
    if (pairspin_ != AlphaBeta) // AlphaAlpha or BetaBeta -> this is a pure triplet
      return Cm_ij_ij;
    if (P==Q) // AlphaBeta C_ii_ii -> this is a pure singlet
      return Cp_ij_ij;
    if(P==O) // AlphaBeta C_ij_ij
      return(0.5*(Cp_ij_ij+Cm_ij_ij));
    else // AlphaBeta C_ij_ji
      return(0.5*(Cp_ij_ij-Cm_ij_ij));
  }
  
}
