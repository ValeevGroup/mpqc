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

#include <vector>
#include <iostream>

#include <util/state/state.h>
#include <util/keyval/keyval.h>
#include <math/scmat/matrix.h>
#include <math/scmat/vector3.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/gaussshell.h>
#include <math/mmisc/grid.h>

#ifdef MPQC_NEW_FEATURES
#include "mpqc/range.hpp"
#endif

using boost::property_tree::ptree;

namespace sc {

class BasisFileSet;
class Integral;

class CartesianIter;
class SphericalTransformIter;

/// @addtogroup ChemistryBasisGaussian
/// @{

/** The GaussianBasisSet class is used describe a basis set composed of atomic gaussian orbitals.
 *  Inputs for common basis sets are included in the MPQC distribution.  They have been obtained
 *  from the EMSL Basis Set Database and translated into the MPQC format.  The citation for this
 *  database is below.  The technical citation for each basis set is listed in the individual
 *  basis set data files, in MPQC's <tt>lib/basis</tt> directory.
 *
 *  Following is a table with available basis sets listing the supported
 *  elements for each basis and the number of basis functions for H, \f$n_0\f$,
 *  first row, \f$n_1\f$, and second row, \f$n_2\f$, atoms.  Basis sets with
 *  non-alpha-numerical characters in their name must be given in quotes.
 *
 *  <table>
 *  <tr><td>Basis Set<td>Elements<td>\f$n_0\f$<td>\f$n_1\f$<td>\f$n_2\f$
 *  <tr><td><tt>STO-2G</tt><td>H-Ca<td>1<td>5<td>9
 *  <tr><td><tt>STO-3G</tt><td>H-Kr<td>1<td>5<td>9
 *  <tr><td><tt>STO-3G*</tt><td>H-Ar<td>1<td>5<td>14
 *  <tr><td><tt>STO-6G</tt><td>H-Kr<td>1<td>5<td>9
 *  <tr><td><tt>MINI (Huzinaga)</tt><td>H-Ca<td>1<td>5<td>9
 *  <tr><td><tt>MINI (Scaled)</tt><td>H-Ca<td>1<td>5<td>9
 *  <tr><td><tt>MIDI (Huzinaga)</tt><td>H-Na, Al-K<td>2<td>9<td>13
 *  <tr><td><tt>DZ (Dunning)</tt><td>H, Li, B-Ne, Al-Cl<td>2<td>10<td>18
 *  <tr><td><tt>DZP (Dunning)</tt><td>H, Li, B-Ne, Al-Cl<td>5<td>16<td>24
 *  <tr><td><tt>DZP + Diffuse (Dunning)</tt><td>H, B-Ne<td>6<td>19<td>
 *  <tr><td><tt>3-21G</tt><td>H-Kr<td>2<td>9<td>13
 *  <tr><td><tt>3-21G*</tt><td>H-Ar<td>2<td>9<td>19
 *  <tr><td><tt>3-21++G</tt><td>H-Ar<td>3<td>13<td>17
 *  <tr><td><tt>3-21++G*</tt><td>H-Ar<td>3<td>13<td>23
 *  <tr><td><tt>4-31G</tt><td>H-Ne, P-Cl<td>2<td>9<td>13
 *  <tr><td><tt>6-31G</tt><td>H-Zn<td>2<td>9<td>13
 *  <tr><td><tt>6-31G*</tt><td>H-Zn<td>2<td>15<td>19
 *  <tr><td><tt>6-31G**</tt><td>H-Zn<td>5<td>15<td>19
 *  <tr><td><tt>6-31+G*</tt><td>H-Ar<td>2<td>19<td>23
 *  <tr><td><tt>6-31++G</tt><td>H-Ca<td>3<td>13<td>17
 *  <tr><td><tt>6-31++G*</tt><td>H-Ar<td>3<td>19<td>23
 *  <tr><td><tt>6-31++G**</tt><td>H-Ar<td>6<td>19<td>23
 *  <tr><td><tt>6-311G</tt><td>H-Ca, Ga-Kr<td>3<td>13<td>21
 *  <tr><td><tt>6-311G*</tt><td>H-Ca, Ga-Kr<td>3<td>18<td>26
 *  <tr><td><tt>6-311G**</tt><td>H-Ca, Ga-Kr<td>6<td>18<td>26
 *  <tr><td><tt>6-311G(2df,2pd)</tt><td>H-Ne, K, Ca<td>14<td>30<td>
 *  <tr><td><tt>6-311++G**</tt><td>H-Ne<td>7<td>22<td>
 *  <tr><td><tt>6-311++G(2d,2p)</tt><td>H-Ca<td>10<td>27<td>35
 *  <tr><td><tt>6-311++G(3df,3pd)</tt><td>H-Ar<td>18<td>39<td>47
 *  <tr><td><tt>cc-pVDZ</tt><td>H-Ar, Ca, Ga-Kr<td>5<td>14<td>18
 *  <tr><td><tt>cc-pVTZ</tt><td>H-Ar, Ca, Ga-Kr<td>14<td>30<td>34
 *  <tr><td><tt>cc-pVQZ</tt><td>H-Ar, Ca, Ga-Kr<td>30<td>55<td>59
 *  <tr><td><tt>cc-pV5Z</tt><td>H-Ar, Ca, Ga-Kr<td>55<td>91<td>95
 *  <tr><td><tt>cc-pV6Z</tt><td>H, He, B-Ne, Al-Ar<td>91<td>140<td>144
 *  <tr><td><tt>aug-cc-pVDZ</tt><td>H, He, B-Ne, Al-Ar, Ga-Kr<td>9<td>23<td>27
 *  <tr><td><tt>aug-cc-pVTZ</tt><td>H, He, B-Ne, Al-Ar, Ga-Kr<td>23<td>46<td>50
 *  <tr><td><tt>aug-cc-pVQZ</tt><td>H, He, B-Ne, Al-Ar, Ga-Kr<td>46<td>80<td>84
 *  <tr><td><tt>aug-cc-pV5Z</tt><td>H, He, B-Ne, Al-Ar, Ga-Kr<td>80<td>127<td>131
 *  <tr><td><tt>aug-cc-pV6Z</tt><td>H, He, B-Ne, Al-Ar<td>127<td>189<td>193
 *  <tr><td><tt>cc-pCVDZ</tt><td>Li, B-Ar<td><td>18<td>27
 *  <tr><td><tt>cc-pCVTZ</tt><td>Li, B-Ar<td><td>43<td>59
 *  <tr><td><tt>cc-pCVQZ</tt><td>Li, B-Ar<td><td>84<td>109
 *  <tr><td><tt>cc-pCV5Z</tt><td>B-Ne<td><td>145<td>
 *  <tr><td><tt>aug-cc-pCVDZ</tt><td>B-F, Al-Ar<td><td>27<td>36
 *  <tr><td><tt>aug-cc-pCVTZ</tt><td>B-Ne, Al-Ar<td><td>59<td>75
 *  <tr><td><tt>aug-cc-pCVQZ</tt><td>B-Ne, Al-Ar<td><td>109<td>134
 *  <tr><td><tt>aug-cc-pCV5Z</tt><td>B-F<td><td>181<td>
 *  <tr><td><tt>NASA Ames ANO</tt><td>H, B-Ne, Al, P, Ti, Fe, Ni<td>30<td>55<td>59
 *  <tr><td><tt>pc-0</tt><td>H, C-F, Si-Cl<td>2<td>9<td>13
 *  <tr><td><tt>pc-1</tt><td>H, C-F, Si-Cl<td>5<td>14<td>18
 *  <tr><td><tt>pc-2</tt><td>H, C-F, Si-Cl<td>14<td>30<td>34
 *  <tr><td><tt>pc-3</tt><td>H, C-F, Si-Cl<td>34<td>64<td>64
 *  <tr><td><tt>pc-4</tt><td>H, C-F, Si-Cl<td>63<td>109<td>105
 *  <tr><td><tt>pc-0-aug</tt><td>H, C-F, Si-Cl<td>3<td>13<td>17
 *  <tr><td><tt>pc-1-aug</tt><td>H, C-F, Si-Cl<td>9<td>23<td>27
 *  <tr><td><tt>pc-2-aug</tt><td>H, C-F, Si-Cl<td>23<td>46<td>50
 *  <tr><td><tt>pc-3-aug</tt><td>H, C-F, Si-Cl<td>50<td>89<td>89
 *  <tr><td><tt>pc-4-aug</tt><td>H, C-F, Si-Cl<td>88<td>145<td>141
 *  </table>
 *
 *  All basis sets
 *  were obtained from the Extensible Computational Chemistry
 *  Environment Basis Set Database, Version 12/03/03, as developed and
 *  distributed by the Molecular Science Computing Facility, Environmental and
 *  Molecular Sciences Laboratory which is part of the Pacific Northwest
 *  Laboratory, P.O. Box 999, Richland, Washington 99352, USA, and funded by
 *  the U.S. Department of Energy. The Pacific Northwest Laboratory is a
 *  multi-program laboratory operated by Battelle Memorial Institute for the
 *  U.S. Department of Energy under contract DE-AC06-76RLO 1830. Contact David
 *  Feller or Karen Schuchardt for further information.
*/
class GaussianBasisSet: virtual public SavableState, virtual public DescribedXMLWritable
{
  public:

    /// Shell is a GaussianShell that is part of GaussianBasisSet, i.e. has a center on which it's centered
    class Shell : public GaussianShell {
      public:
        Shell(const GaussianBasisSet* basis, unsigned int center, const GaussianShell& shell);
        Shell(const Shell& other) :
          GaussianShell(static_cast<const GaussianShell&>(other)),
          basis_(other.basis_),
          center_(other.center_) {}
        ~Shell() {}

        Shell& operator=(const Shell& other) {
          static_cast<GaussianShell&>(*this) = static_cast<const GaussianShell&>(other);
          basis_ = other.basis_;
          center_ = other.center_;
          return *this;
        }

        /// Returns true if this and the argument are equivalent.
        bool equiv(const Shell& s) const;

        const GaussianBasisSet* basis() const { return basis_; }
        unsigned int center() const { return center_; }

        //const double[3]& r() const {}


        virtual ptree& write_xml(ptree& parent, const XMLWriter& writer);

      private:
        friend class GaussianBasisSet;

        static ClassDesc class_desc_;

        const GaussianBasisSet* basis_;
        unsigned int center_;
    };


  private:
    friend class UnionBasisSet;

    ///
    // primary data
    ///
    std::string name_;    // non-empty if keyword "name" was provided
    std::string label_;   // same as name_ if name_ non-empty, else something else
  protected: Ref<Molecule> molecule_;
  private: std::vector<Shell> shells_;

    // computes secondary data from the primary data
    void init2(int skip_ghosts=0,bool include_q=0);

    ///
    // secondary data (expensive or convenient to have precomputed)
    ///
    std::vector<int> shell_to_function_;
    std::vector<int> function_to_shell_;
    std::vector<int> shell_to_primitive_;
    std::vector<int> center_to_shell_;
    std::vector<unsigned int> center_to_nshell_;
    std::vector<unsigned int> center_to_nfunction_;
    unsigned int nbasis_;
    unsigned int nprim_;
    bool has_pure_;
    // old-school SCMatrix stuff, may go away
    Ref<SCMatrixKit> matrixkit_;
    Ref<SCMatrixKit> so_matrixkit_;
    RefSCDimension basisdim_;

    // Counts shells in this basis for this chemical element
    int count_shells_(Ref<KeyVal>& keyval, const char* elemname, const char* sbasisname, BasisFileSet& bases,
		      int havepure, int pure, bool missing_ok);
    // constructs shells on @c atom
    void get_shells_(unsigned int atom, Ref<KeyVal>& keyval, const char* elemname, const char* sbasisname, BasisFileSet& bases,
		     int havepure, int pure, bool missing_ok);
    // Counts shells in an even-tempered primitive basis
    int count_even_temp_shells_(Ref<KeyVal>& keyval, const char* elemname, const char* sbasisname,
                                int havepure, int pure);
    // Constructs an even-tempered primitive basis
    void get_even_temp_shells_(unsigned int atom, Ref<KeyVal>& keyval, const char* elemname, const char* sbasisname,
                               int havepure, int pure);
    // Constructs basis set specified as an array of shells
    void recursively_get_shell(unsigned int atom,Ref<KeyVal>&,
                               const char*,const char*,BasisFileSet&,
                               int,int,int,bool missing_ok);

    // finishes up the KeyVal constructor
    void init(Ref<Molecule>&,Ref<KeyVal>&,
              BasisFileSet&,
              int have_userkeyval,
              int pure);

    /**
     * verifies that symmetry-equivalent atoms have the same basis set
     * this may help to squash hard-to-detect bugs due to the user providing an erroneous basis keyword
     * will throw InputError if the basis set does not conform with the point group
     */
    void validate_point_group() const;

  protected:
    /* This CTOR leaves it up to the derived class to completely
        initialize the basis set. */
    GaussianBasisSet();
    virtual void set_matrixkit(const Ref<SCMatrixKit>&);

    /** Initializes everything. To be used by derived classes. */
    void init(std::string name,
              std::string label,
              const Ref<Molecule> &molecule,
              const std::vector<Shell>& shell);

  public:
    /** This holds scratch data needed to compute basis function values. */
    class ValueData {
      protected:
        CartesianIter **civec_;
        SphericalTransformIter **sivec_;
        int maxam_;
      public:
        ValueData(const Ref<GaussianBasisSet> &, const Ref<Integral> &);
        ~ValueData();
        CartesianIter **civec() { return civec_; }
        SphericalTransformIter **sivec() { return sivec_; }
    };

    /// @name Constructors
    ///@{

    /** The KeyVal constructor.

        <dl>

        <dt><tt>molecule</tt><dd> The gives a Molecule object.  The is no
        default.

        <dt><tt>puream</tt><dd> If this boolean parameter is true then 5D,
        7F, etc. will be used.  Otherwise all cartesian functions will be
        used.  The default depends on the particular basis set.

        <dt><tt>name</tt><dd> This is a string giving the name of the basis
        set.  The above table of basis sets gives some of the recognized
        basis set names.  It may be necessary to put the name in double
        quotes. There is no default.

        <dt><tt>basis</tt><dd> If the <tt>element</tt> vector is given,
        then this vector specifies the names of basis sets
        for each element. If the <tt>element</tt> vector is not given,
        this vector specifies basis set name for each atom in the molecule
        (note that the same basis name must be specified
        for each set of atoms related by symmetry).
        If this keyword omitted, the basis
        set specified in <tt>name</tt> will be used for all atoms.

        <dt><tt>element</tt><dd> This is a vector of elements.  If it is
        given then it must have the same number of entries as the basis
        vector.

        <dt><tt>basisdir</tt><dd> A string giving a directory where basis
        set data files are to be sought.  See the text below for a complete
        description of what directories are consulted.

        <dt><tt>basisfiles</tt><dd> Each keyword in this vector of files is
        appended to the directory specified with basisdir and basis set
        data is read from them.

        </dl>

        Several files in various directories are checked for basis set
        data.  First, basis sets can be given by the user in the basis
        section at the top level of the main input file.  Next, if a path
        is given with the basisdir keyword, then all of the files given
        with the basisfiles keyword are read in after appending their names
        to the value of basisdir.  Basis sets can be given in these files
        in the basis section at the top level as well.  If the named basis
        set still cannot be found, then GaussianBasisSet will try convert
        the basis set name to a file name (see the note below for the rules
        of this conversion) and check first in the directory
        given by basisdir.  Next it checks for the environment variable
        MPQC_DATA_PATH.  If it is set it will look for the basis file in
        $MPQC_DATA_PATH/basis.  Otherwise it will look in the source code
        distribution in the directory SC/lib/basis.  If the executable has
        changed machines or the source code has be moved, then it may be
        necessary to copy the library files to your machine and set the
        MPQC_DATA_PATH environmental variable.

        <b>Note</b>: translation of a basis name to a file name will convert
        upper-case letters(A-Z) to the lower-case letters,
        characters ',' and ' ' (whitespace) to '_', character '+' to 'P',
        character '*' to 'S', character '(' to 'L', and character ')' to 'R'.

        The basis set itself is also given in the ParsedKeyVal format. There are two
        recognized formats for basis sets:
        <dl>

        <dt>array of shells<dd> One must specify the keyword :basis: followed by the
        lowercase atomic name followed by : followed by the basis set name
        (which may need to be placed inside double quotes). The value for the keyword
        is an array of shells. Each shell
        reads the following keywords:

        <dl>

        <dt><tt>type</tt><dd> This is a vector that describes each
        component of this shell.  For each element the following two
        keywords are read:

        <dl>

          <dt><tt>am</tt><dd> The angular momentum of the component.  This
          can be given as the letter designation, s, p, d, etc.  There is
          no default.

          <dt><tt>puream</tt><dd> If this boolean parameter is true then
          5D, 7F, etc. shells are used.  The default is false.  This
          parameter can be overridden in the GaussianBasisSet
          specification.

        </dl>

        <dt><tt>exp</tt><dd> This is a vector giving the exponents of the
        primitive Gaussian functions.

        <dt><tt>coef</tt><dd> This is a matrix giving the coeffients of the
        primitive Gaussian functions.  The first index gives the component
        number of the shell and the second gives the primitive number.

        </dl>

        <dt><dd>An example might be easier to understand.  This is a basis set
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

        <dt>basis set of even-tempered primitive Gaussians<dd>
        Such basis set format is given as a group of keywords. The name of the group is :basis: followed by the
        lowercase atomic name followed by : followed by the basis set name
        (which may need to be placed inside double quotes).
        The group of keywords must contain vectors <tt>am</tt> and <tt>nprim</tt>,
        which specify the angular momentum and the number of shells in each
        block of even-tempered primitives. In addition, one must provide any
        two of the following vectors:

        <dl>
          <dt><tt>first_exp</tt><dd> The exponent of the "tightest" primitive Gaussian in the block.

          <dt><tt>last_exp</tt><dd> The exponent of the most "diffuse" primitive Gaussian in the block.

          <dt><tt>exp_ratio</tt><dd> The ratio of exponents of consecutive primitive Gaussians
          in the block.

        </dl>

        <dt><dd>Note that the dimensions of each vector must be the same.

        Here's an example of a basis set composed of 2 blocks of even-tempered s-functions
        and 1 block of even-tempered p-functions.

        <pre>
        basis: (
         neon: "20s5s13p":(

           am = [ 0 0 1 ]
           nprim = [ 20 5 13 ]
           first_exp = [ 1000.0 0.1  70.0 ]
           last_exp =  [    1.0 0.01  0.1 ]

         )
        )
        </pre>

        </dl>

        */
    GaussianBasisSet(const Ref<KeyVal>&);

    /**
     * Constructs GaussianBasisSet from a Molecule and a vector of GaussianShells
     * @param molecule molecule whose centers will serve as origins for @c shells
     * @param shells vector of GaussianShell objects
     * @param shell_to_center maps shell -> center; shells on the same center <b>MUST</b> appear in succession
     * @param name
     * @param label
     */
    GaussianBasisSet(const Ref<Molecule>& molecule,
                     const std::vector<GaussianShell>& shells,
                     const std::vector<unsigned int>& shell_to_center,
                     std::string name = std::string(),
                     std::string label = std::string());

    GaussianBasisSet(const GaussianBasisSet&);

    /** This will produce a GaussianBasisSet object composed of a a single "unit" basis function,
     * i.e. that is one everywhere. This can be used with integral evaluators to compute certain classes
     * of integrals, with limitations. */
    static Ref<GaussianBasisSet> unit();

    ///@}

    /// @name Serialization
    ///@{

    /// restores this from @c si
    GaussianBasisSet(StateIn& si);
    /// saves this to @c so
    void save_data_state(StateOut& so);

    virtual ptree& write_xml(ptree& parent, const XMLWriter& writer);

    ///@}

    virtual ~GaussianBasisSet();

    /** Uses UnionBasisSet to construct the sum of this and B.  */
    Ref<GaussianBasisSet> operator+(const Ref<GaussianBasisSet>& B);
    GaussianBasisSet& operator=(const GaussianBasisSet& A);

    /// Return the name of the basis set (is nonnull only if keyword "name" was provided)
    const std::string& name() const { return name_; }
    /** Return the label of the basis set. label() returns the same string as name() if
        keyword "name" was provided, otherwise it is a unique descriptive string which
        can be arbitrarily long. */
    const std::string& label() const { if (!name_.empty()) return name_; else return label_; }

    /// Return the Molecule object.
    Ref<Molecule> molecule() const { return molecule_; }
    /// Returns the SCMatrixKit that is to be used for AO bases.
    Ref<SCMatrixKit> matrixkit() { return matrixkit_; }
    /// Returns the SCMatrixKit that is to be used for SO bases.
    Ref<SCMatrixKit> so_matrixkit() { return so_matrixkit_; }
    /// Returns the SCDimension object for the dimension.
    RefSCDimension basisdim() { return basisdim_; }

    /// Return the number of centers.
    unsigned int ncenter() const;
    /// Return the number of shells.
    unsigned int nshell() const { return shells_.size(); }
    /// Return the number of shells on the given center.
    int nshell_on_center(int icenter) const;
    /** Return an overall shell number, given a center and the shell number
        on that center. \sa shell_to_center() */
    int shell_on_center(int icenter, int shell) const;
    /// Return the center on which the given shell is located. This is a non-decreasing function, i.e.
    /// if s1 > s2 then shell_to_center(s1) >= shell_to_center(s2)
    int shell_to_center(int ishell) const { return shells_[ishell].center(); }
    /// Return the number of basis functions.
    unsigned int nbasis() const { return nbasis_; }
    /// Return the number of basis functions on the given center.
    unsigned int nbasis_on_center(int icenter) const;
    /// Return the number of primitive Gaussians (sum of number of primitives for each Shell)
    unsigned int nprimitive() const { return nprim_; }
    /// Return true if basis contains solid harmonics Gaussians
    bool has_pure() const { return has_pure_; }

    /// Return the maximum number of functions that any shell has.
    unsigned int max_nfunction_in_shell() const;
    /** Return the maximum number of Cartesian functions that any shell has.
        The optional argument is an angular momentum increment. */
    unsigned int max_ncartesian_in_shell(int aminc=0) const;
    /// Return the maximum number of primitive Gaussian that any shell has.
    unsigned int max_nprimitive_in_shell() const;
    /// Return the highest angular momentum in any shell.
    unsigned int max_angular_momentum() const;
    /// Return the maximum number of Gaussians in a contraction in any shell.
    unsigned int max_ncontraction() const;
    /** Return the maximum angular momentum found in the given contraction
        number for any shell.  */
    unsigned int max_am_for_contraction(int con) const;
    /// Return the maximum number of Cartesian functions in any shell.
    unsigned int max_cartesian() const;

    /// Return the number of the first function in the given shell.
    int shell_to_function(int i) const { return shell_to_function_[i]; }
    /// Return the shell to which the given function belongs.
    int function_to_shell(int i) const;

#ifdef MPQC_NEW_FEATURES
    /// @return range object for the basis functions in shell s
    mpqc::range range(int s) const;
#endif

    /// Return a reference to Shell number i.
    const Shell& operator()(int i) const { return shells_[i]; }
    /// Return a reference to Shell number i.
    Shell& operator()(int i) { return shells_[i]; }
    /// Return a reference to Shell number i.
    const Shell& operator[](int i) const { return shells_[i]; }
    /// Return a reference to Shell number i.
    Shell& operator[](int i) { return shells_[i]; }
    /// Return a reference to Shell number i.
    const Shell& shell(int i) const { return shells_[i]; }
    /// Return a reference to Shell number i.
    Shell& shell(int i) { return shells_[i]; }

    /// Return a reference to Shell number ishell on center icenter.
    const Shell& operator()(int icenter,int ishell) const;
    /// Return a reference to Shell number ishell on center icenter.
    Shell& operator()(int icenter,int ishell);
    /// Return a reference to Shell number j on center i.
    const Shell& shell(int i,int j) const { return operator()(i,j); }
    /// Return a reference to Shell number j on center i.
    Shell& shell(int i,int j) { return operator()(i,j); }

    /// Return the absolute index of shell S located at center C in this basis. If the shell is not found, returns -1
    int find(int C, const GaussianShell& S) const;

    /** The location of center icenter.  The xyz argument is 0 for x, 1 for
        y, and 2 for z. */
    double r(int icenter,int xyz) const;

    /** Compute the values for this basis set at position r.  The
        basis_values argument must be vector of length nbasis. */
    int values(const SCVector3& r, ValueData *, double* basis_values) const;
    /** Like values(...), but computes gradients of the basis function
        values, too.  The g_values argument must be vector of length
        3*nbasis.  The data will be written in the order bf1_x, bf1_y,
        bf1_z, ... */
    int grad_values(const SCVector3& r, ValueData *,
                    double*g_values,double* basis_values=0) const;
    /** Like values(...), but computes first and second derivatives of the
        basis function values, too.  h_values must be vector of length
        6*nbasis.  The data will be written in the order bf1_xx, bf1_yx,
        bf1_yy, bf1_zx, bf1_zy, bf1_zz, ... */
    int hessian_values(const SCVector3& r, ValueData *, double *h_values,
                       double*g_values=0,double* basis_values=0) const;
    /** Compute the values for the given shell functions at position r.
        See the other values(...) members for more information.  */
    int shell_values(const SCVector3& r, int sh,
                     ValueData *, double* basis_values) const;
    /** Like values(...), but computes gradients of the shell function
        values, too.  See the other grad_values(...)
        members for more information.  */
    int grad_shell_values(const SCVector3& r, int sh,
                          ValueData *,
                          double*g_values, double* basis_values=0) const;
    /** Like values(...), but computes first and second derivatives of the
        shell function values, too.  See the other hessian_values(...)
        members for more information. */
    int hessian_shell_values(const SCVector3& r, int sh,
                       ValueData *, double *h_values,
                       double*g_values=0,double* basis_values=0) const;

    /// Returns true if this and the argument are equivalent.
    int equiv(const Ref<GaussianBasisSet> &b);

    /// Print a brief description of the basis set.
    void print_brief(std::ostream& =ExEnv::out0()) const;
    /// Print a detailed description of the basis set.
    void print(std::ostream& =ExEnv::out0()) const;
};

/// Nonmember operator+ is more convenient to use than the member operator+
Ref<GaussianBasisSet>
operator+(const Ref<GaussianBasisSet>& A, const Ref<GaussianBasisSet>& B);

/** Find A in bs on center icenter and return its index */
int
ishell_on_center(int icenter, const Ref<GaussianBasisSet>& bs,
	             const GaussianShell& A);

/// Nonmember operator+ is more convenient to use than the member operator+
Ref<GaussianBasisSet>
operator+(const Ref<GaussianBasisSet>& A, const Ref<GaussianBasisSet>& B);

/** computes a map from basis functions in A to the equivalent basis functions in B. A must be contained in B,
    else will throw */
std::vector<unsigned int> operator<<(const GaussianBasisSet& B, const GaussianBasisSet& A);
/// same as operator<<, except A does not have to be contained in B, map[a] = -1 if function a is not in B
std::vector<int> map(const GaussianBasisSet& B, const GaussianBasisSet& A);

/// A heavy-duty map from one GaussianBasisSet to another GaussianBasisSet. Can only be constructed
/// if the original basis set is a subset of the target basis set.
class GaussianBasisSetMap : public RefCount {
  public:
    /// will throw if Source is not contained in Target
    GaussianBasisSetMap(const Ref<GaussianBasisSet>& Source, const Ref<GaussianBasisSet>& Target);
    ~GaussianBasisSetMap();

    const Ref<GaussianBasisSet>& source() const { return source_; }
    const Ref<GaussianBasisSet>& target() const { return target_; }

    /// maps shell s in Source to its location in Target
    int map_shell(int s) const;
    /// maps function f in Source to its location in Target
    int map_function(int f) const;
    /** it is useful to organize sets of basis functions that are adjacent in Source
        and in Target. These sets will be termed fblocks. Thus each fblock is has size of at least 1 and there
        is at least 1 fblock. This returns the number of fblocks. */
    int nfblock() const;
    /// returns the Source index of the first basis function in fblock b
    int fblock_to_function(int b) const;
    /// the number of basis functions in fblock b
    int fblock_size(int b) const;

  private:
    Ref<GaussianBasisSet> source_;
    Ref<GaussianBasisSet> target_;

    /////////
    // maps
    /////////
    /// maps shells in Source to their location in Target
    std::vector<int> smap_;
    /// maps functions in Source to their location in Target
    std::vector<int> fmap_;
    /// return the first bf in fblock (using Source indexing)
    std::vector<int> fblock_to_function_;
    /// return the size of fblock
    std::vector<int> fblock_size_;

};

class WriteBasisGrid : public WriteVectorGrid {
  private:
    struct BasisFunctionMap : public DimensionMap {
      std::vector<int> map;
      int operator()(int o) const { return map[o]; }
      std::size_t ndim() const { return map.size(); }
    };
    Ref<GaussianBasisSet> basis_;
    Ref<Integral> integral_;

  protected:
    BasisFunctionMap bfmap_;

    void label(char* buffer);
    Ref<Molecule> get_molecule() { return basis_->molecule(); }
    void calculate_values(const std::vector<SCVector3>& Points, std::vector<double>& Values);
    const DimensionMap& dimension_map() const { return bfmap_; }
    std::size_t ndim() const { return bfmap_.ndim(); }
    void initialize();

    static Ref<KeyVal> process_keyval_for_base_class(const Ref<KeyVal>& kv);

  public:

    WriteBasisGrid(const Ref<KeyVal>& keyval);
    WriteBasisGrid(
        const Ref<GaussianBasisSet>& basis,
        const Ref<Grid>& grid,
        std::string gridformat,
        std::string gridfile
    );

    virtual ~WriteBasisGrid();

    virtual ptree& write_xml(ptree& parent, const XMLWriter& writer);

};

} // end namespace sc
/// @}
// end of addtogroup ChemistryBasisGaussian

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
