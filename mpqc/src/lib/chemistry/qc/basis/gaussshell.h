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

#include <iostream>
#include <util/state/state.h>
#include <math/scmat/vector3.h>
#include <util/keyval/keyval.h>

namespace sc {

class CartesianIter;
class SphericalTransformIter;
class Integral;

/// A Gaussian orbital shell.
class GaussianShell: public SavableState
{
  public:
    enum PrimitiveType { Normalized, Unnormalized };
    enum GaussianType { Cartesian, Pure };
  private:
    int nprim;
    int ncon;
    int* l;
    int* puream;
    double* exp;
    double** coef;  // contraction coefficients for unnormalized primitives

    // computed data:
    int nfunc;
    int min_am_;
    int max_am_;
    int ncart_;
    int has_pure_;
    int* contr_to_func_;
    int* func_to_contr_;
    void init_computed_data();

    double shell_normalization(int);
    void convert_coef();
    void normalize_shell();
    PrimitiveType keyval_init(const Ref<KeyVal>&,int,int);
    int test_monobound(double &r, double &bound) const;
  public:

    static const char* amtypes;
    static const char* AMTYPES;

    /** A GaussianShell constructor.
        Users of GaussianShell must pass pointers to newed memory that is kept
        by GaussianShell and deleted by the destructor.
        The arguments for the following ctor are:
        <ul>
          <li> ncn is the number of contracted functions
             (1 except for SP and gen. con.)
          <li> nprm is the number of primitives
          <li> e gives the exponents (length nprm)
          <li> am gives the angular momentum (length ncn)
          <li> pure is 1 for pure am and 0 for cartesian (length ncn)
          <li> c are the contraction coefficients (length ncn by nprm)
          <li> pt describes whether the primitive functions are to be
            considered normalized or unnormalized.  This effects whether
            or not c is manipulated to give the correct normalization.
          <li> If do_normalize_shell is true (the default), then the
            shell normalization constants will be folded into the
            coefficients.
        </ul>
    */
    GaussianShell(
                  int ncn,
                  int nprm,
                  double* e,
                  int* am,
                  int* pure,
                  double** c,
                  PrimitiveType pt = GaussianShell::Normalized,
                  bool do_normalize_shell = true);
    /** A GaussianShell constructor.  In this ctor pure is either
     GaussianShell::Cartesian or Gaussian::Pure and all of the contracted
     functions are treated in that way. (The user doesn\'t need to compute
     generate a int*pure vector in this case.) */
    GaussianShell(
                  int ncn,
                  int nprm,
                  double* e,
                  int* am,
                  GaussianType pure,
                  double** c,
                  PrimitiveType pt = GaussianShell::Normalized);
    /// Construct a GaussianShell from KeyVal input.
    GaussianShell(const Ref<KeyVal>&);
    /// Restore a GaussianShell from a StateIn object.
    GaussianShell(StateIn&);
    /** Construct a GaussianShell from KeyVal input.  If pure
        is nonzero Cartesian functions will be used, otherwise,
        solid harmonics will be used. */
    GaussianShell(const Ref<KeyVal>&,int pure);
    ~GaussianShell();
    void save_data_state(StateOut&);
    /// The number of primitive Gaussian shells.
    int nprimitive() const { return nprim; }
    /// The number of contractions formed from the primitives.
    int ncontraction() const { return ncon; }
    /// The number of basis functions.
    int nfunction() const { return nfunc; }
    /// The maximum angular momentum in the shell.
    int max_angular_momentum() const { return max_am_; }
    /// The minimum angular momentum in the shell.
    int min_angular_momentum() const { return min_am_; }
    /// The maximum number of Cartesian functions in any contraction.
    int max_cartesian() const;
    /// The angular momentum of the given contraction.
    int am(int con) const { return l[con]; }
    /// The maximum angular momentum of any contraction.
    int max_am() const { return max_am_; }
    /// The minimum angular momentum of any contraction.
    int min_am() const { return min_am_; }
    /// The character symbol for the angular momentum of the given contraction.
    char amchar(int con) const { return amtypes[l[con]]; }
    /// The number of basis functions coming from the given contraction.
    int nfunction(int con) const;
    /// The total number of functions if this shell was Cartesian.
    int ncartesian() const { return ncart_; }
    /** The total number of Cartesian functions if this shift is applied to
        all of the angular momentums. */
    int ncartesian_with_aminc(int aminc) const;
    /// The number of Cartesian functions for the given contraction.
    int ncartesian(int con) const { return ((l[con]+2)*(l[con]+1))>>1; }
    /// Returns nonzero if contraction con is Cartesian.
    int is_cartesian(int con) const { return !puream[con]; }
    /// Returns nonzero if contraction con is solid harmonics.
    int is_pure(int con) const { return puream[con]; }
    /// Returns nonzero if any contraction is solid harmonics.
    int has_pure() const { return has_pure_; }
    /// Returns the number of the first function in the given contraction
    int contraction_to_function(int c) const { return contr_to_func_[c]; }
    /// Returns the contraction to which this function belongs
    int function_to_contraction(int f) const { return func_to_contr_[f]; }
    /// Returns the contraction coef for unnormalized primitives.
    double coefficient_unnorm(int con,int prim) const {return coef[con][prim];}
    /// Returns the contraction coef for normalized primitives.
    double coefficient_norm(int con,int prim) const;
    /// Returns the exponent of the given primitive.
    double exponent(int iprim) const { return exp[iprim]; }

    /** Compute the values for this shell at position r.  The
        basis_values argument must be vector of length nfunction(). */
    int values(CartesianIter **, SphericalTransformIter **,
               const SCVector3& r, double* basis_values);
    /** Like values(...), but computes gradients of the basis function
        values, too. */
    int grad_values(CartesianIter **, SphericalTransformIter **,
                    const SCVector3& R,
                    double* g_values,
                    double* basis_values=0) const;
    /** Like values(...), but computes first and second derivatives of the
        basis function values, too. */
    int hessian_values(CartesianIter **, SphericalTransformIter **,
                       const SCVector3& R,
                       double* h_values, double* g_values=0,
                       double* basis_values=0) const;

    /** Returns the intra-generalized-contraction overlap
        matrix element <con func1|con func2> within an arbitrary
        constant for the shell. */
    double relative_overlap(const Ref<Integral>&,
                            int con, int func1, int func2) const;
    /** Returns the intra-generalized-contraction overlap matrix element
        <con func1|con func2> within an arbitrary constant for the shell.
        func1 and func2 are determined according to the axis exponents, a1,
        b1, c1, a2, b2, and c2. */
    double relative_overlap(int con,
                            int a1, int b1, int c1,
                            int a2, int b2, int c2) const;

    /// Returns true if this and the argument are equivalent.
    int equiv(const GaussianShell* s) const;
    /// Returns true if this and the argument are equivalent.
    int equiv(const GaussianShell& s) const
    {
      return equiv(&s);
    }


    /** Returns a radius.  All functions in the shell are below
        threshold outside this radius. */
    double extent(double threshold) const;

    /** Returns a bound for the basis function.  This bound
        is defined so that it is positive and monotonically
        decreasing as a function of r. */
    double monobound(double r) const;

    void print(std::ostream& =ExEnv::out0()) const;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
