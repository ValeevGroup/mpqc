//
// inttypedescr.h
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
#pragma interface
#endif

#include <util/ref/ref.h>

#ifndef _chemistry_qc_basis_inttypedescr_h
#define _chemistry_qc_basis_inttypedescr_h

namespace sc {
  
  /** For a kind of integrals (e.g. <ij|1/r12|kl> )
      IntegralTypeDescr gives basic properties with respect to permutations, etc.
  */
  class IntegralTypeDescr : public RefCount {
  public:
    IntegralTypeDescr() {}
    virtual ~IntegralTypeDescr() {}

    virtual unsigned int num_particles() const =0;
    /// number of functions for particle i
    virtual unsigned int num_functions(unsigned int i) const =0;
    /// Symmetry with respect to permutation of functions for particle i
    virtual int perm_symm(unsigned int i) const =0;
    /// Symmetry with respect to permutation of particles i and j
    virtual int perm_symm(unsigned int i, unsigned int j) const =0;
  };

  /** Two-body integrals with n functions per particle and given symmetry properties */
  class TwoBodyIntTypeDescr : public IntegralTypeDescr {
  public:
    TwoBodyIntTypeDescr(unsigned int n, int perm_p1, int perm_p2, int perm_p12);

    /// Implementation of IntegralTypeDescr::num_particles()
    unsigned int num_particles() const;
    /// Implementation of IntegralTypeDescr::num_functions()
    unsigned int num_functions(unsigned int i) const;
    /// Implementation of IntegralTypeDescr::perm_symm()
    int perm_symm(unsigned int i) const;
    /// Implementation of IntegralTypeDescr::perm_symm()
    int perm_symm(unsigned int i, unsigned int j) const;

  private:
    unsigned int n_;
    int perm_p1_;
    int perm_p2_;
    int perm_p12_;

  };
  
}

#endif

