//
// gaussbas.h
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

#ifndef _chemistry_qc_basis_gaussbas_h
#define _chemistry_qc_basis_gaussbas_h

#ifdef __GNUC__
#pragma interface
#endif

#include <iostream.h>

#include <util/state/state.h>
#include <util/state/array.h>
#include <math/scmat/matrix.h>
#include <math/scmat/vector3.h>
#include <chemistry/molecule/molecule.h>

class GaussianShell;
class RefKeyVal;
class BasisFileSet;

SavableState_REF_fwddec(GaussianBasisSet);
SavableState_REF_fwddec(Integral)

class CartesianIter;
class SphericalTransformIter;

/** The GaussianBasisSet class is used describe a basis set composed of
atomic gaussian orbitals.

Following is a table with available basis sets listing the supported elements
for each basis and the number of basis functions for H, $n_0$, first row,
$n_1$, and second row, $n_2$, atoms:

\begin{tabular}{lcrrr}
  \multicolumn{1}{c}{Basis Set}&
   \multicolumn{1}{c}{Elements}&
   \multicolumn{1}{c}{$n_0$}&
   \multicolumn{1}{c}{$n_1$}&
   \multicolumn{1}{c}{$n_2$} \\
\verb*|STO-2G| & H-Ca & 1 & 5 & 9 \\
\verb*|STO-3G| & H-Kr & 1 & 5 & 9 \\
\verb*|STO-3G*| & H-Ar & 1 & 5 & 15 \\
\verb*|STO-6G| & H-Kr & 1 & 5 & 9 \\
\verb*|MINI (Huzinaga)| & H-Ca & 1 & 5 & 9 \\
\verb*|MINI (Scaled)| & H-Ca & 1 & 5 & 9 \\
\verb*|MIDI (Huzinaga)| & H-Na & 2 & 9 &  \\
\verb*|DZ (Dunning)| & H, Li, B-Ne, Al-Cl & 2 & 10 & 18 \\
\verb*|DZP (Dunning)| & H, Li, B-Ne, Al-Cl & 5 & 16 & 24 \\
\verb*|DZP + Diffuse (Dunning)| & H, B-Ne & 6 & 19 &  \\
\verb*|3-21G| & H-Kr & 2 & 9 & 13 \\
\verb*|3-21G*| & H-Ar & 2 & 9 & 19 \\
\verb*|3-21++G| & H-Ar & 3 & 13 & 17 \\
\verb*|3-21++G*| & H-Ar & 3 & 13 & 23 \\
\verb*|4-31G| & H-Ne, P-Cl & 2 & 9 & 13 \\
\verb*|4-31G*| & H-Ne, P-Cl & 2 & 15 & 19 \\
\verb*|4-31G**| & H-Ne, P-Cl & 5 & 15 & 19 \\
\verb*|6-31G| & H-Ar & 2 & 9 & 13 \\
\verb*|6-31G*| & H-Ar & 2 & 15 & 19 \\
\verb*|6-31G**| & H-Ar & 5 & 15 & 19 \\
\verb*|6-31+G*| & H-Ar & 2 & 19 & 23 \\
\verb*|6-31++G| & H-Ar & 3 & 13 & 17 \\
\verb*|6-31++G*| & H-Ar & 3 & 19 & 23 \\
\verb*|6-31++G**| & H-Ar & 6 & 19 & 23 \\
\verb*|6-311G| & H-Ar, Ga-Kr & 3 & 13 & 21 \\
\verb*|6-311G*| & H-Ar, Ga-Kr & 3 & 19 & 27 \\
\verb*|6-311G**| & H-Ar, Ga-Kr & 6 & 19 & 27 \\
\verb*|6-311G(2df,2pd)| & H-Ne & 15 & 35 &  \\
\verb*|6-311++G**| & H-Ne & 7 & 23 &  \\
\verb*|6-311++G(2d,2p)| & H-Ne & 10 & 29 &  \\
\verb*|6-311++G(3df,3pd)| & H-Ar & 19 & 45 & 53 \\
\verb*|cc-pVDZ| & H, He, B-Ne, Al-Ar & 5 & 14 & 18 \\
\verb*|cc-pVTZ| & H, He, B-Ne, Al-Ar & 14 & 30 & 34 \\
\verb*|cc-pVQZ| & H, He, B-Ne, Al-Ar & 30 & 55 & 59 \\
\verb*|cc-pV5Z| & H-Ne, Al-Ar & 55 & 91 & 95 \\
\verb*|aug-cc-pVDZ| & H, He, B-Ne, Al-Ar & 9 & 23 & 27 \\
\verb*|aug-cc-pVTZ| & H, He, B-Ne, Al-Ar & 23 & 46 & 50 \\
\verb*|aug-cc-pVQZ| & H, He, B-Ne, Al-Ar & 46 & 80 & 84 \\
\verb*|aug-cc-pV5Z| & H, He, B-Ne, Al-Ar & 80 & 127 & 131 \\
\verb*|cc-pCVDZ| & B-Ne &  & 18 &  \\
\verb*|cc-pCVTZ| & B-Ne &  & 43 &  \\
\verb*|cc-pCVQZ| & B-Ne &  & 84 &  \\
\verb*|cc-pCV5Z| & B-Ne &  & 145 &  \\
\verb*|aug-cc-pCVDZ| & B-F &  & 27 &  \\
\verb*|aug-cc-pCVTZ| & B-Ne &  & 59 &  \\
\verb*|aug-cc-pCVQZ| & B-Ne &  & 109 &  \\
\verb*|aug-cc-pCV5Z| & B-F &  & 181 &  \\
\verb*|NASA Ames ANO| & H, B-Ne, Al, P, Ti, Fe, Ni & 30 & 55 & 59 \\
\end{tabular}

*/
class GaussianBasisSet: public SavableState
{
#   define CLASSNAME GaussianBasisSet
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    char* name_;
    GaussianShell** shell_;
    Arrayint shell_to_function_;
    Arrayint function_to_shell_;

    RefMolecule molecule_;

    RefSCMatrixKit matrixkit_;
    RefSCMatrixKit so_matrixkit_;
    RefSCDimension basisdim_;

    // these are needed for the routines to compute basis set values
    // they must be initialized with set_integral()
    CartesianIter **civec_;
    SphericalTransformIter **sivec_;

    int ncenter_;
    SSBArrayint shell_to_center_;
    SSBArrayint center_to_shell_;
    SSBArrayint center_to_nshell_;
    SSBArrayint center_to_nbasis_;

    int nshell_;
    int nbasis_;
    int nprim_;

    void recursively_get_shell(int&,RefKeyVal&,
                               const char*,const char*,BasisFileSet&,
                               int,int,int);

    void init(RefMolecule&,RefKeyVal&,
              BasisFileSet&,
              int have_userkeyval,
              int pure);
    void init2(int skip_ghosts=0);
    
  protected:
    GaussianBasisSet(const GaussianBasisSet&);
    virtual void set_matrixkit(const RefSCMatrixKit&);
    
  public:
    /** @memo The KeyVal constructor.

        \begin{description}

        \item[molecule] The gives a Molecule object.  The is no default.

        \item[puream] If this boolean parameter is true then 5D, 7F,
        etc. will be used.  Otherwise all cartesian functions will be used.
        The default depends on the particular basis set.

        \item[name] This is a string giving the name of the basis set.  The
        above table of basis sets gives some of the recognized basis set
        names.  It may be necessary to put the name in double quotes. There
        is no default.

        \item[basis] This is a vector of basis set names that can give a
        different basis set to each atom in the molecule.  If the element
        vector is given, then it gives different basis sets to different
        elements.  The default is to give every atom the basis set
        specified in name.

        \item[element] This is a vector of elements.  If it is given then
        it must have the same number of entries as the basis vector.

        \item[basisdir] A string giving a directory where basis set data
        files are to be sought.  See the text below for a complete
        description of what directors are consulted.

        \item[basisfiles] Each keyword in this vector of files is appended
        to the directory specified with basisdir and basis set data is read
        from them.

        \item[matrixkit] Specifies a SCMatrixKit object.  It is usually not
        necessary to give this keyword, as the default action should get
        the correct SCMatrixKit.

        \end{description}

        Several files in various directories are checked for basis set
        data.  First, basis sets can be given by the user in the basis
        section at the top level of the main input file.  Next, if a path
        is given with the basisdir keyword, then all of the files given
        with the basisfiles keyword are read in after appending their names
        to the value of basisdir.  Basis sets can be given in these files
        in the basis section at the top level as well.  If the named basis
        set still cannot be found, then GaussianBasisSet will try convert
        the basis set name to a file name and check first in the directory
        given by basisdir.  Next it checks for the environment variable
        SCLIBDIR.  If it is set it will look for the basis file in
        \$SCLIBDIR/basis.  Otherwise it will look in the source code
        distribution in the directory SC/lib/basis.  If the executable has
        changed machines or the source code has be moved, then it may be
        necessary to copy the library files to your machine and set the
        SCLIBDIR environmental variable.

        The basis set itself is also given in the ParsedKeyVal format.  It
        is a vector of shells with the keyword :basis: followed by the
        lowercase atomic name followed by : followed by the basis set name
        (which may need to be placed inside double quotes).  Each shell
        reads the following keywords:

        \begin{description}

        \item[type] This is a vector that describes each
        component of this shell.  For each element the following two
        keywords are read:

        \begin{description}

          \item[am] The angular momentum of the component.  This can be
          given as the letter designation, s, p, d, etc.  There is no
          default.

          \item[puream] If this boolean parameter is true then 5D, 7F,
          etc. shells are used.  The default is false.  This parameter can
          be overridden in the GaussianBasisSet specification.

        \end{description}

        \item[exp] This is a vector giving the exponents of the primitive
        Gaussian functions.

        \item[coef] This is a matrix giving the coeffients of the primitive
        Gaussian functions.  The first index gives the component number of
        the shell and the second gives the primitive number.

        \end{description}

        An example might be easier to understand.  This is a basis set
        specificition for STO-2G carbon:

        <pre>
        basis: (
         carbon: "STO-2G": [
          (type: [(am = s)]
           {      exp      coef:0 } = {
              27.38503303 0.43012850
               4.87452205 0.67891353
           })
          (type: [(am = p) (am = s)]
           {     exp      coef:1     coef:0 } = {
               1.13674819 0.04947177 0.51154071
               0.28830936 0.96378241 0.61281990
           })
         ]
        )
        </pre>
     */
    GaussianBasisSet(const RefKeyVal&);
    GaussianBasisSet(StateIn&);
    virtual ~GaussianBasisSet();

    void save_data_state(StateOut&);

    /// Return the name of the basis set.
    const char* name() const { return name_; }

    /// Return the Molecule object.
    RefMolecule molecule() const { return molecule_; }
    /// Returns the SCMatrixKit that is to be used for AO bases.
    RefSCMatrixKit matrixkit() { return matrixkit_; }
    /// Returns the SCMatrixKit that is to be used for SO bases.
    RefSCMatrixKit so_matrixkit() { return so_matrixkit_; }
    /// Returns the SCDimension object for the dimension.
    RefSCDimension basisdim() { return basisdim_; }

    /// Return the number of centers.
    int ncenter() const;
    /// Return the number of shells.
    int nshell() const { return nshell_; }
    /// Return the number of shells on the given center.
    int nshell_on_center(int icenter) const;
    /** Return an overall shell number, given a center and the shell number
        on that center. */
    int shell_on_center(int icenter, int shell) const;
    /// Return the center on which the given shell is located.
    int shell_to_center(int ishell) const { return shell_to_center_(ishell); }
    /// Return the number of basis functions.
    int nbasis() const { return nbasis_; }
    /// Return the number of basis functions on the given center.
    int nbasis_on_center(int icenter) const;
    /// Return the number of primitive Gaussians.
    int nprimitive() const { return nprim_; }

    /// Return the maximum number of functions that any shell has.
    int max_nfunction_in_shell() const;
    /** Return the maximum number of Cartesian functions that any shell has.
        The optional argument is an angular momentum increment. */
    int max_ncartesian_in_shell(int aminc=0) const;
    /// Return the highest angular momentum in any shell.
    int max_angular_momentum() const;
    /// Return the maximum number of Gaussians in a contraction in any shell.
    int max_ncontraction() const;
    /** Return the maximum angular momentum found in the given contraction
        number for any shell.  */
    int max_am_for_contraction(int con) const;
    /// Return the maximum number of Cartesian functions in any shell.
    int max_cartesian() const;

    /// Return the number of the first function in the given shell.
    int shell_to_function(int i) const { return shell_to_function_(i); }
    /// Return the shell to which the given function belongs.
    int function_to_shell(int i) const;

    /// Return a reference to GaussianShell number i.
    const GaussianShell& operator()(int i) const { return *shell_[i]; }
    /// Return a reference to GaussianShell number i.
    GaussianShell& operator()(int i) { return *shell_[i]; }
    /// Return a reference to GaussianShell number i.
    const GaussianShell& operator[](int i) const { return *shell_[i]; }
    /// Return a reference to GaussianShell number i.
    GaussianShell& operator[](int i) { return *shell_[i]; }
    /// Return a reference to GaussianShell number i.
    const GaussianShell& shell(int i) const { return *shell_[i]; }
    /// Return a reference to GaussianShell number i.
    GaussianShell& shell(int i) { return *shell_[i]; }

    /// Return a reference to GaussianShell number ishell on center icenter.
    const GaussianShell& operator()(int icenter,int ishell) const;
    /// Return a reference to GaussianShell number ishell on center icenter.
    GaussianShell& operator()(int icenter,int ishell);
    /// Return a reference to GaussianShell number j on center i.
    const GaussianShell& shell(int i,int j) const { return operator()(i,j); }
    /// Return a reference to GaussianShell number j on center i.
    GaussianShell& shell(int i,int j) { return operator()(i,j); }

    /** The location of center icenter.  The xyz argument is 0 for x, 1 for
        y, and 2 for z. */
    double r(int icenter,int xyz) const;
    
    /** Compute the values for this basis set at position r.  The
        basis_values argument must be vector of length nbasis. */
    int values(const SCVector3& r, double* basis_values) const;
    /** Like values(...), but computes gradients of the basis function
        values, too.  The g_values argument must be vector of length
        3*nbasis.  The data will be written in the order bf1_x, bf1_y,
        bf1_z, ... */
    int grad_values(const SCVector3& r,
                    double*g_values,double* basis_values=0) const;
    /** Like values(...), but computes first and second derivatives of the
        basis function values, too.  h_values must be vector of length
        6*nbasis.  The data will be written in the order bf1_xx, bf1_yx,
        bf1_yy, bf1_zx, bf1_zy, bf1_zz, ... */
    int hessian_values(const SCVector3& r, double *h_values,
                       double*g_values=0,double* basis_values=0) const;
    /** Compute the values for the given shell functions at position r.
        See the other values(...) members for more information.  */
    int shell_values(const SCVector3& r, int sh, double* basis_values) const;
    /** Like values(...), but computes gradients of the shell function
        values, too.  See the other grad_values(...)
        members for more information.  */
    int grad_shell_values(const SCVector3& r, int sh,
                          double*g_values, double* basis_values=0) const;
    /** Like values(...), but computes first and second derivatives of the
        shell function values, too.  See the other hessian_values(...)
        members for more information. */
    int hessian_shell_values(const SCVector3& r, int sh, double *h_values,
                       double*g_values=0,double* basis_values=0) const;
    /** This must be called before the values, grid_values, and
        hessian_values members to initialize iterators that know the basis
        function order. */
    void set_integral(const RefIntegral&);

    /// Returns true if this and the argument are equivalent.
    int equiv(const RefGaussianBasisSet &b);

    /// Print a brief description of the basis set.
    void print_brief(ostream& =cout) const;
    /// Print a detailed description of the basis set.
    void print(ostream& =cout) const;
};

SavableState_REF_dec(GaussianBasisSet);

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
