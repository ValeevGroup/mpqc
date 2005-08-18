//
// utils.h
//
// Copyright (C) 2005 Edward Valeev
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

#ifdef __GNUG__
#pragma interface
#endif

#include <vector>
#include <math/scmat/matrix.h>

#ifndef _chemistry_qc_mbptr12_utils_h
#define _chemistry_qc_mbptr12_utils_h

namespace sc {
  
  class MOIndexSpace;
  
  /** Antisymmetrizes 4-index quantity <ij|A|kl> -> <ij|A|kl> - <ij|A|lk>
      and saves to Aanti. Row dimension has to be an integer multiple of
      bra->rank()*bra->rank(). Same for ket.
    */
  void antisymmetrize(RefSCMatrix& Aanti, const RefSCMatrix& A,
                      const Ref<MOIndexSpace>& bra,
                      const Ref<MOIndexSpace>& ket,
                      bool accumulate = false);
  
  /** Converts RefDiagSCMatrix to std::vector<double>
  */
  std::vector<double> convert(const RefDiagSCMatrix& A);
  
  /// print out the Fortran-style matrix
  void print_f77_mat(const std::string& comment,
                     const double* A,
                     unsigned int nrow,
                     unsigned int ncol,
                     bool transpose = false);
}

#endif

