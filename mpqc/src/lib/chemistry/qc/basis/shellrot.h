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

class Integral;
class ShellRotation {
  private:
    int n_;
    int am_;
    double **r;
    
    void done();

  public:
    void init(int a, SymmetryOperation&, const RefIntegral&);
    void init_pure(int a, SymmetryOperation&, const RefIntegral&);
    
    ShellRotation(int n);
    ShellRotation(const ShellRotation&);
    ShellRotation(int a, SymmetryOperation&, const RefIntegral&, int pure =0);
    virtual ~ShellRotation();

    ShellRotation& operator=(const ShellRotation&);
    
    int am() const { return am_; }
    int dim() const { return n_; }
    
    double& operator()(int i, int j) { return r[i][j]; }
    double* operator[](int i) { return r[i]; }
    
    ShellRotation operate(const ShellRotation&) const;
    ShellRotation sim_transform(const ShellRotation&) const;
    
    double trace() const;
    
    void print() const;
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
