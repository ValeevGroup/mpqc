//
// inttypedescr.cc
//
// Copyright (C) 2007 Edward Valeev
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
#pragma implementation
#endif

#include <util/class/scexception.h>
#include <chemistry/qc/basis/inttypedescr.h>

using namespace sc;

////

TwoBodyIntTypeDescr::TwoBodyIntTypeDescr(unsigned int n, int perm_p1, int perm_p2, int perm_p12) :
  n_(n), perm_p1_(perm_p1), perm_p2_(perm_p2), perm_p12_(perm_p12)
{
}

unsigned int TwoBodyIntTypeDescr::num_particles() const { return 2; }
unsigned int TwoBodyIntTypeDescr::num_functions(unsigned int i) const { return n_; }
int TwoBodyIntTypeDescr::perm_symm(unsigned int i) const {
  if (i == 1) return perm_p1_;
  if (i == 2) return perm_p2_;
  throw ProgrammingError("IntBodyIntTypeDescr::perm_symm(i) -- i must be 1 or 2",__FILE__,__LINE__);
}
int TwoBodyIntTypeDescr::perm_symm(unsigned int i, unsigned int j) const {
  if ( (i == 1 && j == 2) ||
       (i == 2 && j == 1) ) return perm_p12_;
  throw ProgrammingError("IntBodyIntTypeDescr::perm_symm(i,j) -- i,j must be 1,2",__FILE__,__LINE__);
}
