//
// shellrot.h
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

#ifndef _chemistry_qc_basis_shellrot_h
#define _chemistry_qc_basis_shellrot_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/symmetry/pointgrp.h>

namespace sc {

class Integral;
/** Compute the transformation matrices that maps a set of Cartesian
    functions to another set of Cartesian functions in a rotated
    coordinate system. */
class ShellRotation {
  private:
    int n_;
    int am_;
    double **r;
    
    void done();

  public:
    /** Initialize the ShellRotation for Cartesian functions, given the
        angular momentum, a symmetry operation, and an Integral object. */
    void init(int a, SymmetryOperation&, const Ref<Integral>&);
    /** Initialize the ShellRotation for solid harmonic functions, given
        the angular momentum, a symmetry operation, and an Integral
        object. */
    void init_pure(int a, SymmetryOperation&, const Ref<Integral>&);
    
    /// Initialize this ShellRotation to hold a n by n transformation.
    ShellRotation(int n);
    /// Initialize this from another ShellRotation.
    ShellRotation(const ShellRotation&);
    /// Initialize using init(...) or, if pure is nonzero, init_pure(...).
    ShellRotation(int a, SymmetryOperation&, const Ref<Integral>&, int pure =0);
    virtual ~ShellRotation();

    /// Assign this to another shell rotation.
    ShellRotation& operator=(const ShellRotation&);
    
    /// Return the angular momentum.
    int am() const { return am_; }
    /// Return the number of functions in a shell.
    int dim() const { return n_; }
    
    /// Return an element of the transform matrix.
    double& operator()(int i, int j) { return r[i][j]; }
    /// Return a row of the transform matrix.
    double* operator[](int i) { return r[i]; }
    
    /// Returns the result of rot*this.
    ShellRotation operate(const ShellRotation&rot) const;
    /// Returns the result of rot*this*transpose(rot).
    ShellRotation transform(const ShellRotation&rot) const;
    
    /// Return the trace of the transformation.
    double trace() const;
    
    /// Print the object to ExEnv::out0().
    void print() const;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
