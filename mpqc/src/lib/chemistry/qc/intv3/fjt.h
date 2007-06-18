//
// fjt.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
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

#ifndef _chemistry_qc_intv3_fjt_h
#define _chemistry_qc_intv3_fjt_h

#include <chemistry/qc/basis/fjt.h>

namespace sc {

/// Computes F_j(T) using 6-th order Taylor interpolation
class FJT: public Fjt {
  private:
    double **gtable;

    int maxj;
    double *denomarray;
    double wval_infinity;
    int itable_infinity;

    double *int_fjttable;

    int ngtable() const { return maxj + 7; }
  public:
    FJT(int n);
    ~FJT();
    /// implementation of Fjt::values()
    double *values(int J, double T);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
