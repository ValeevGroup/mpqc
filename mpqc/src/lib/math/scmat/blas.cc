//
// blas.cc
//
// Copyright (C) 2007 Edward Valeev
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

#include <cmath>
#include <stdexcept>
#include <util/misc/formio.h>
#include <math/scmat/blas.h>

using namespace std;
using namespace sc;

namespace sc {

  void C_DGEMM(char transa, char transb, int m, int n, int k, double alpha,
               const double *A, int nca, const double *B, int ncb, double beta,
               double *C, int ncc)
  {

    /* the only strange thing we need to do is reverse everything
       since the stride runs differently in C vs. Fortran
     */
    
    /* also, do nothing if a dimension is 0 */
    if (m == 0 || n == 0 || k == 0) return;

    F77_DGEMM(&transb,&transa,&n,&m,&k,&alpha,B,&ncb,A,&nca,&beta,C,&ncc);

  }
  
};

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
