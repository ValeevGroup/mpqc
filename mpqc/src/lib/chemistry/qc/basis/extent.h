//
// extent.h
//
// Copyright (C) 2000 Sandia National Labs
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
// Maintainer: SNL
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

#ifndef _chemistry_qc_basis_extent_h
#define _chemistry_qc_basis_extent_h

#ifdef __GNUC__
#pragma interface
#endif

#include <float.h>
#include <chemistry/qc/basis/basis.h>

struct ExtentData {
    int shell;
    double bound;
    ExtentData() {}
    ExtentData(int s, double b): shell(s), bound(b) {}
};

ARRAY_dec(ExtentData);

class ShellExtent: public VRefCount {
    double lower_[3];
    double resolution_;
    int n_[3];
    ArrayExtentData *contributing_shells_;
    ArrayExtentData null_;

    ArrayExtentData &data(int *b);
    double distance(double loc, int axis, int origin, int point);
    ArrayExtentData &data(int x, int y, int z);
  public:
    ShellExtent();
    ~ShellExtent();
    void init(const RefGaussianBasisSet&,
              double resolution = 1.0, double tolerance = DBL_EPSILON);
    /** Returns the shells that are nonzero at coordinates x, y, z.
        The shells numbers are in ascending order. */
    const ArrayExtentData &contributing_shells(int x, int y, int z)
        { return data(x,y,z); }
    const ArrayExtentData &contributing_shells(double x, double y, double z);
    void print(ostream &o = ExEnv::out());
    const int *n() const { return n_; }
    int n(int ixyz) const { return n_[ixyz]; }
    double lower(int ixyz) const { return lower_[ixyz]; }
    double resolution() const { return resolution_; }
};
REF_dec(ShellExtent);

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
