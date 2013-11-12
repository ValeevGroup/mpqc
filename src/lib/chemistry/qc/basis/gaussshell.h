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

#include <iostream>
#include <util/state/state.h>
#include <math/scmat/vector3.h>
#include <util/keyval/keyval.h>
#include <util/misc/xml.h>

using boost::property_tree::ptree;

namespace sc {

class CartesianIter;
class SphericalTransformIter;
class Integral;

/// A shell of Gaussian functions. A shell is a set of functions with same quantum numbers, contraction coefficients,
/// and exponents, and located on the common origin. GaussianShell does include the origin information.
/// @sa GaussianBasisSet::Shell
class GaussianShell: public DescribedXMLWritable
{
  public:
    enum PrimitiveType { Normalized, Unnormalized };
    enum GaussianType { Cartesian, Pure };
    /// @return numerical epsilon value; exponents/coefficients that differ by less than this are the same,
    /// similarly coefficients less than that in absolute magntide are zero (zero exponents are OK!!!)
    static double epsilon() { return 1.0e-13; }
  private:
    // primary data
    std::vector<unsigned int> l;  //!< "angular momenta" of contracted shells, vector of ncontr elements
    std::vector<bool> puream;     //!< pure(true) or cart(false) angular character, vector of ncontr elements
    std::vector<double> exp;      //!< exponents, vector of nprim elements
    std::vector<double> coef_blk; //!< contraction coefficients for unnormalized primitives, vector of ncontr times nprim elements

    // computes secondary data; may also chomp coefs and exponents using epsilon(), hence may modify exp and coef_blk
    void init_computed_data();
    static void chomp(std::vector<double>& exp, std::vector<double>& coef_blk,
                      unsigned int ncontr, double epsilon = GaussianShell::epsilon());

    // secondary data:
    std::vector<double*> coef;  //!< contraction coefficients accessable as vector of pointers
    int nfunc;                  //!< number of basis functions
    int min_am_;                //!< the lowest "angular momentum" of any contraction in the shell
    int max_am_;                //!< the highest "angular momentum" of any contraction in the shell
    int ncart_;                 //!< number of Cartesian basis functions, same as nfunc if all elements of @c puream are false
    int has_pure_;
    int has_cartesian_;
    std::vector<int> contr_to_func_;
    std::vector<int> func_to_contr_;

    double shell_normalization(int);
    void convert_coef();
    void normalize_shell();
    PrimitiveType keyval_init(const Ref<KeyVal>&,int,int);
    int test_monobound(double &r, double &bound) const;

    friend void ToStateOut(const GaussianShell &s, StateOut &so, int &count);
    friend void FromStateIn(GaussianShell &s, StateIn &si, int &count);

  public:

    // Default constructor needed to work with sc::SavableState
    GaussianShell(){}

    static const char* amtypes;
    static const char* AMTYPES;

    /**
     * Constructs GaussianShell programmatically. The constructed shells has @c ncontr contractions
     * composed of @c nprim primitive Gaussians.
     * @param am angular momenta of each contracted shell, std::vector of @c ncontr nonnegative integers
     * @param puream spherical/cartesian flag, std::vector of @c ncontr boolean values
     * @param exps exponents of each primitive, std::vector of @c nprim nonnegative double-precision numbers
     * @param contr_coefs contraction coefficients for each contraction, std::vector of @c ncontr times @c nprim
     *                    double-precision numbers. @c contr_coefs[c*ncontr+p] is the contraction coefficient
     *                    of primitive @c p in contraction @c c
     * @param pt indicates whether the input contraction coefficients are in terms of normalized primitive functions (i.e.
     *           do not include the normalization factors in them).
     * @param normalize_shell if true, the coefficients will be scaled to normalize each contraction to unity.
     */
    GaussianShell(const std::vector<unsigned int>& am,
                  const std::vector<bool>& puream,
                  const std::vector<double>& exps,
                  const std::vector<double>& contr_coefs,
                  PrimitiveType pt = GaussianShell::Normalized,
                  bool normalize_shell = true);

    /** Construct a GaussianShell object from KeyVal input.
     *
     * @param kv the KeyVal object
     * @param pure this is an optional parameter that can be used to override programmatically
     *             the "pure" keyword value in @c kv. If @c pure=0 Cartesian GaussianShell will be
     *             constructed; if @c pure=1 solid(spherical) harmonics GaussianShell will be
     *             constructed; any other value (or omitting this value) will use the @c kv object.
     */
    GaussianShell(const Ref<KeyVal>& kv, int pure=-1);
    /// Copy constructor (deep :-)
    GaussianShell(const GaussianShell& other);

    /** @deprecated
     *  A GaussianShell constructor.
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
          <li> c are the contraction coefficients (C-style array of arrays!! not single block ncn by nprm)
          <li> pt describes whether the primitive functions are to be
            considered normalized or unnormalized.  This effects whether
            or not c is manipulated to give the correct normalization.
          <li> If do_normalize_shell is true (the default), then the
            shell normalization constants will be folded into the
            coefficients.
        </ul>
    */
    DEPRECATED GaussianShell(
                  int ncn,
                  int nprm,
                  double* e,
                  int* am,
                  int* pure,
                  double** c,
                  PrimitiveType pt = GaussianShell::Normalized,
                  bool do_normalize_shell = true);
    /** @deprecated A GaussianShell constructor.  In this ctor pure is either
     GaussianShell::Cartesian or Gaussian::Pure and all of the contracted
     functions are treated in that way. (The user doesn\'t need to compute
     generate a int*pure vector in this case.) */
    DEPRECATED GaussianShell(
                  int ncn,
                  int nprm,
                  double* e,
                  int* am,
                  GaussianType pure,
                  double** c,
                  PrimitiveType pt = GaussianShell::Normalized);

    /**
     *
     * @return the "unit" Gaussian shell, zero angular momentum, and zero exponents and unit contraction coefficient
     */
    static GaussianShell unit();

    ~GaussianShell();
    /// The number of primitive Gaussian shells.
    unsigned int nprimitive() const { return exp.size(); }
    /// The number of contractions formed from the primitives.
    unsigned int ncontraction() const { return l.size(); }
    /// The number of basis functions.
    unsigned int nfunction() const { return nfunc; }
    /// The maximum angular momentum in the shell.
    int max_angular_momentum() const { return max_am_; }
    /// The minimum angular momentum in the shell.
    int min_angular_momentum() const { return min_am_; }
    /// The maximum number of Cartesian functions in any contraction.
    int max_cartesian() const;
    /// The angular momenta of contractions
    const std::vector<unsigned int>& am() const { return l; }
    /// The angular momentum of the given contraction.
    unsigned int am(int con) const { return l[con]; }
    /// The maximum angular momentum of any contraction.
    unsigned int max_am() const { return max_am_; }
    /// The minimum angular momentum of any contraction.
    unsigned int min_am() const { return min_am_; }
    /// The character symbol for the angular momentum of the given contraction.
    char amchar(int con) const { return amtypes[l[con]]; }
    /// The number of basis functions coming from the given contraction.
    unsigned int nfunction(int con) const;
    /// The total number of functions if this shell was Cartesian.
    unsigned int ncartesian() const { return ncart_; }
    /** The total number of Cartesian functions if this shift is applied to
        all of the angular momentums. */
    unsigned int ncartesian_with_aminc(int aminc) const;
    /// The number of Cartesian functions for the given contraction.
    unsigned int ncartesian(int con) const { return ((l[con]+2)*(l[con]+1))>>1; }
    /// Returns nonzero if contraction con is Cartesian.
    bool is_cartesian(int con) const { return !puream[con]; }
    /// Returns nonzero if any contraction is Cartesian.
    bool has_cartesian() const { return has_cartesian_; }
    /// Returns true if contraction con is solid harmonics.
    bool is_pure(int con) const { return puream[con]; }
    /// Vector of booleans that indicate whether each contraction is solid harmonics.
    const std::vector<bool>& is_pure() const { return puream; }
    /// Returns true if any contraction is solid harmonics.
    bool has_pure() const { return has_pure_; }
    /// Returns the number of the first function in the given contraction
    int contraction_to_function(int c) const { return contr_to_func_[c]; }
    /// Returns the contraction to which this function belongs
    int function_to_contraction(int f) const { return func_to_contr_[f]; }
    /// Returns the contraction coef for unnormalized primitives.
    double coefficient_unnorm(int con,int prim) const {return coef[con][prim];}
    /// Returns the contraction coef for normalized primitives.
    double coefficient_norm(int con,int prim) const;
    /// returns coefficients for unnormalization primitives, in block form
    const std::vector<double>& coefficient_unnorm_block() const { return coef_blk; }
    /// Returns the exponents of the given primitive.
    double exponent(int iprim) const { return exp[iprim]; }
    /// Returns the exponents.
    const std::vector<double>& exponents() const { return this->exp; }

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
    bool equiv(const GaussianShell& s) const;


    /** Returns a radius.  All functions in the shell are below
        threshold outside this radius. */
    double extent(double threshold) const;

    /** Returns a bound for the basis function.  This bound
        is defined so that it is positive and monotonically
        decreasing as a function of r. */
    double monobound(double r) const;

    void print(std::ostream& =ExEnv::out0()) const;

    virtual ptree& write_xml(
        ptree& parent, const XMLWriter& writer
    );
};

  /** constructs a new GaussianShell from @c shell by applying Filter @c filter
   * @tparam Filter boolear functor that returns true if contraction is to be kept,
   *                must provide bool Filter(const GaussianShell& shell, unsigned int contr)
   * @param shell
   * @param f Filter object
   * @return
   */
  template <typename Filter>
  GaussianShell filter(const GaussianShell& shell,
                       Filter f) {

    std::vector<unsigned int> am;
    std::vector<bool> pure;

    for(unsigned int c=0; c<shell.ncontraction(); ++c) {
      if (f(shell, c)) {
        am.push_back(shell.am(c));
        pure.push_back(shell.is_pure(c));
      }
    }

    // filter the coefficients also
    const size_t nprim = shell.nprimitive();
    std::vector<double> coefs(am.size() * shell.nprimitive());
    for(size_t c=0, cp=0; c<shell.ncontraction(); ++c) {
      if (f(shell, c)) {
        std::copy(shell.coefficient_unnorm_block().begin() + c*nprim,
                  shell.coefficient_unnorm_block().begin() + (c+1)*nprim,
                  coefs.begin() + cp);
        cp += nprim;
      }
    }

    return GaussianShell(am, pure, shell.exponents(), coefs, GaussianShell::Unnormalized);
  }

/// writes GaussianShell to sc::StateOut
void ToStateOut(const GaussianShell &s, StateOut &so, int &count);

/// reads GaussianShell from sc::StateIn
void FromStateIn(GaussianShell &s, StateIn &si, int &count);

} // namespace sc

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
