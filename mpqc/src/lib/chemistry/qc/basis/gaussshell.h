//
// gaussshell.h
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

#ifndef _chemistry_qc_basis_gaussshell_h
#define _chemistry_qc_basis_gaussshell_h

#ifdef __GNUC__
#pragma interface
#endif

#include <iostream.h>
#include <util/state/state.h>
#include <math/scmat/vector3.h>

class RefKeyVal;
SavableState_REF_fwddec(Integral)

class GaussianShell: public SavableState
{
#   define CLASSNAME GaussianShell
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    enum PrimitiveType { Normalized, Unnormalized };
    enum GaussianType { Cartesian, Pure };
  private:
    int nprim;
    int ncon;
    int nfunc;
    int* l;
    int* puream;
    double* exp;
    double** coef;  // contraction coefficients for unnormalized primitives

    double shell_normalization(int);
    void convert_coef();
    void normalize_shell();
    void compute_nfunc();
    PrimitiveType keyval_init(const RefKeyVal&,int,int);
    static const char* amtypes;
    static const char* AMTYPES;
  public:
    // Users of GaussianShell must pass pointers to newed memory that is kept
    // by GaussianShell and deleted by the destructor.
    // The arguments for the following ctor are:
    //   ncn is the number of contracted functions
    //     (1 except for SP and gen. con.)
    //   nprm is the number of primitives
    //   e gives the exponents (length nprm)
    //   am gives the angular momentum (length ncn)
    //   pure is 1 for pure am and 0 for cartesian (length ncn)
    //   c are the contraction coefficients (length ncn by nprm)
    //   pt describes whether the primitive functions are to be considered
    //     normalized or unnormalized.  this effects whether or not c is
    //     manipulated to give the correct normalization.
    GaussianShell(
                  int ncn,
                  int nprm,
                  double* e,
                  int* am,
                  int* pure,
                  double** c,
                  PrimitiveType pt = GaussianShell::Normalized);
    // In this ctor pure is either GaussianShell::Cartesian or
    // Gaussian::Pure and all of the contracted functions are
    // treated in that way. (The user doesn\'t need to compute
    // generate a int*pure vector in this case.)
    GaussianShell(
                  int ncn,
                  int nprm,
                  double* e,
                  int* am,
                  GaussianType pure,
                  double** c,
                  PrimitiveType pt = GaussianShell::Normalized);
    GaussianShell(const RefKeyVal&);
    GaussianShell(StateIn&);
    GaussianShell(const RefKeyVal&,int pure);
    ~GaussianShell();
    void save_data_state(StateOut&);
    int nprimitive() const { return nprim; }
    int ncontraction() const { return ncon; }
    int nfunction() const { return nfunc; }
    int max_angular_momentum() const;
    int min_angular_momentum() const;
    int max_cartesian() const;
    int am(int con) const { return l[con]; }
    int max_am() const { return max_angular_momentum(); }
    int min_am() const { return min_angular_momentum(); }
    char amchar(int con) const { return amtypes[l[con]]; }
    int nfunction(int con) const;
    int ncartesian() const;
    // this is given a shift for all of the angular momentums
    int ncartesian_with_aminc(int aminc) const;
    int ncartesian(int con) const;
    int is_cartesian(int con) const { return !puream[con]; }
    int is_pure(int con) const { return puream[con]; }
    int has_pure() const;
    // returns the con coef for unnormalized primitives
    double coefficient_unnorm(int con,int prim) const {return coef[con][prim];}
    // returns the con coef for normalized primitives
    double coefficient_norm(int con,int prim) const;
    double exponent(int iprim) const { return exp[iprim]; }

    // compute the value of this shell at offset r
    int values(const RefIntegral&, const SCVector3& r, double* basis_values);
    int grad_values(const RefIntegral&, const SCVector3& R,
                    double* g_values,
                    double* basis_values=0) const;

    // returns the intra-generalized-contraction overlap
    // matrix element <con func1|con func2> within an arbitrary
    // constant for the shell
    double relative_overlap(const RefIntegral&,
                            int con, int func1, int func2) const;
    double relative_overlap(int con,
                            int a1, int b1, int c1,
                            int a2, int b2, int c2) const;

    void print(ostream& =cout) const;
};

SavableState_REF_dec(GaussianShell);

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
