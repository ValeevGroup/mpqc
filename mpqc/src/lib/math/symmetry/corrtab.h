//
// corrtab.h
//
// Copyright (C) 1997 Limit Point Systems, Inc.
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

#ifdef __GNUC__
#pragma interface
#endif

#ifndef _math_symmetry_corrtab_h
#define _math_symmetry_corrtab_h

#include <iostream.h>

#include <math/symmetry/pointgrp.h>

////////////////////////////////////////////////////////////////////
 
//.  The \clsnm{CorrelationTable} class provides a correlation
//table between two point groups.
class CorrelationTable: public VRefCount {
  private:
    RefPointGroup group_;
    RefPointGroup subgroup_;

    int n_;
    int *ngamma_;
    int **gamma_;
  public:
    //. Create a correlation table for the two groups.
    CorrelationTable(const RefPointGroup& group,
                     const RefPointGroup& subgroup);

    ~CorrelationTable();

    //. Returns the number of irreps in the high order group.
    int n() const { return n_; }
    //. Returns the number of irreps in the low order group that an irrep
    //from the high order group can be reduced to.
    int ngamma(int igamma) const { return ngamma_[igamma]; }
    //. Returns the irreps in the low order group that an irrep from the
    //high order group can be reduced to.
    int gamma(int igamma, int i) const { return gamma_[igamma][i]; }

    void print(ostream &o=cout) const;
};
REF_dec(CorrelationTable);

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
