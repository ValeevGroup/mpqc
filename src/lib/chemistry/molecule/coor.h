//
// coor.h
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

#ifndef _chemistry_molecule_coor_h
#define _chemistry_molecule_coor_h

#include <iostream>
#include <vector>

#include <math/scmat/matrix.h>
#include <math/optimize/transform.h>
#include <chemistry/molecule/molecule.h>

namespace sc {

  /// @addtogroup ChemistryMolecule
  /// @{

/** The IntCoor abstract class describes an internal coordinate of a
molecule. */
class IntCoor: public SavableState {
  protected:
    // conversion factors from radians, bohr to the preferred units
    static double bohr_conv;
    static double radian_conv;
    char *label_;
    double value_;
  public:
    IntCoor(StateIn&);
    IntCoor(const IntCoor&);
    /** This constructor takes a string containing a label for the
        internal coordinate.  The string is copied. */
    IntCoor(const char* label = 0);
    /** The KeyVal constructor.
        <dl>

        <dt><tt>label</tt><dd> A label for the coordinate using only to
        identify the coordinate to the user in printouts.  The default is
        no label.

        <dt><tt>value</tt><dd> A value for the coordinate.  In the way that
        coordinates are usually used, the default is to compute a value
        from the cartesian coordinates in a Molecule object.

        <dt><tt>unit</tt><dd> The unit in which the value is given.  This
        can be bohr, anstrom, radian, and degree.  The default is bohr for
        lengths and radian for angles.

        </dl> */
    IntCoor(const Ref<KeyVal>&);
    
    virtual ~IntCoor();
    void save_data_state(StateOut&);

    /// Returns the string containing the label for the internal coordinate.
    virtual const char* label() const;
    /// Returns the value of the coordinate in atomic units or radians.
    virtual double value() const;
    /// Sets the value of the coordinate in atomic units or radians.
    virtual void set_value(double);
    /// Returns the value of the coordinate in more familiar units.
    virtual double preferred_value() const;
    /// Returns a string representation of the type of coordinate this is.
    virtual const char* ctype() const = 0;
    /// Print information about the coordinate.
    virtual void print(std::ostream & o=ExEnv::out0()) const;
    virtual void print_details(const Ref<Molecule> &, std::ostream& =ExEnv::out0()) const;
    /** Returns the value of the force constant associated with this
        coordinate. */
    virtual double force_constant(Ref<Molecule>&) = 0;
    /// Recalculate the value of the coordinate.
    virtual void update_value(const Ref<Molecule>&) = 0;
    /// Fill in a row the the B matrix.
    virtual void bmat(const Ref<Molecule>&,RefSCVector&bmat,double coef=1.0) = 0;
    /** Test to see if this internal coordinate is equivalent to that one.
        The definition of equivalence is left up to the individual
        coordinates. */
    virtual int equivalent(Ref<IntCoor>&) = 0;
};

/** SumIntCoor is used to construct linear combinations of internal
coordinates.

The following is a sample ParsedKeyVal input for a SumIntCoor object:
<pre>
  sumintcoor\<SumIntCoor>: (
    coor: [
      \<StreSimpleCo>:( atoms = [ 1 2 ] )
      \<StreSimpleCo>:( atoms = [ 2 3 ] )
      ]
    coef = [ 1.0 1.0 ]
    )
</pre>
*/
class SumIntCoor: public IntCoor {
  private:
    std::vector<double> coef_;
    std::vector<Ref<IntCoor> > coor_;
  public:
    SumIntCoor(StateIn&);
    /** This constructor takes a string containing a label for this
        coordinate. */
    SumIntCoor(const char *);
    /** The KeyVal constructor.
        <dl>

        <dt><tt>coor</tt><dd> A vector of IntCoor objects that define the
        summed coordinates.

        <dt><tt>coef</tt><dd> A vector of floating point numbers that gives
        the coefficients of the summed coordinates.

        </dl> */
    SumIntCoor(const Ref<KeyVal>&);

    ~SumIntCoor();
    void save_data_state(StateOut&);

    /// Returns the number of coordinates in this linear combination.
    int n();
    /** Add a coordinate to the linear combination.  coef is the
        coefficient for the added coordinate. */
    void add(Ref<IntCoor>&,double coef);
    /// This function normalizes all the coefficients.
    void normalize();

    // IntCoor overrides
    /// Returns the value of the coordinate in a.u. and radians.
    double preferred_value() const;
    /// Always returns ``SUM''.
    const char* ctype() const;
    /// Print the individual coordinates in the sum with their coefficients.
    void print_details(const Ref<Molecule> &, std::ostream& =ExEnv::out0()) const;
    /// Returns the weighted sum of the individual force constants.
    double force_constant(Ref<Molecule>&);
    /// Recalculate the value of the coordinate.
    void update_value(const Ref<Molecule>&);
    /// Fill in a row the the B matrix.
    void bmat(const Ref<Molecule>&,RefSCVector&bmat,double coef = 1.0);
    /// Always returns 0.
    int equivalent(Ref<IntCoor>&);
};

/** The SetIntCoor class describes a set of internal coordinates.
It can automatically generate these coordinates using a integral coordinate
generator object (see the IntCoorGen class) or the internal
coordinates can be explicity given.

The following is a sample ParsedKeyVal input for
a SetIntCoor object.
<pre>
  setintcoor<SetIntCoor>: [
    \<SumIntCoor>: (
      coor: [
        \<StreSimpleCo>:( atoms = [ 1 2 ] )
        \<StreSimpleCo>:( atoms = [ 2 3 ] )
        ]
      coef = [ 1.0 1.0 ]
      )
    \<BendSimpleCo>:( atoms = [ 1 2 3 ] )
  ]
</pre>
*/
class SetIntCoor: public SavableState {
  private:
    std::vector<Ref<IntCoor> > coor_;
  public:
    SetIntCoor();
    SetIntCoor(StateIn&);
    /** The KeyVal constructor.
        <dl>

        <dt><tt>generator</tt><dd> A IntCoorGen object that will be used to
        generate the internal coordinates.

        <dt><tt>i</tt><dd> A sequence of integer keywords, all \f$i\f$ for
        \f$0 \leq i < n\f$, can be assigned to IntCoor objects.

        </dl> */
    SetIntCoor(const Ref<KeyVal>&);

    virtual ~SetIntCoor();
    void save_data_state(StateOut&);

    /// Adds an internal coordinate to the set.
    void add(const Ref<IntCoor>&);
    /// Adds all the elements of another set to this one.
    void add(const Ref<SetIntCoor>&);
    /// Removes the last coordinate from this set.
    void pop();
    /// Removes all coordinates from the set.
    void clear();
    /// Returns the number of coordinates in the set.
    int n() const;
    /// Returns a reference to the i'th coordinate in the set.
    Ref<IntCoor> coor(int i) const;
    /// Compute the B matrix by finite displacements.
    virtual void fd_bmat(const Ref<Molecule>&,RefSCMatrix&);
    /// Compute the B matrix the old-fashioned way.
    virtual void bmat(const Ref<Molecule>&, RefSCMatrix&);
    /** Create an approximate Hessian for this set of coordinates.  This
        Hessian is a symmetric matrix whose i'th diagonal is the force
        constant for the i'th coordinate in the set. */
    virtual void guess_hessian(Ref<Molecule>&,RefSymmSCMatrix&);
    /// Print the coordinates in the set.
    virtual void print_details(const Ref<Molecule> &,std::ostream& =ExEnv::out0()) const;
    /// Recalculate the values of the internal coordinates in the set.
    virtual void update_values(const Ref<Molecule>&);
    /// Copy the values of the internal coordinates to a vector.
    virtual void values_to_vector(const RefSCVector&);
};


// ////////////////////////////////////////////////////////////////////////

class BitArrayLTri;

/** IntCoorGen generates a set of simple internal coordinates
    for a molecule. */
class IntCoorGen: public SavableState
{
  protected:
    Ref<Molecule> molecule_;
    
    int linear_bends_;
    int linear_lbends_;
    int linear_tors_;
    int linear_stors_;
    int nextra_bonds_;
    int *extra_bonds_;
    double linear_bend_thres_;
    double linear_tors_thres_;
    double radius_scale_factor_;

    void init_constants();

    double cos_ijk(Molecule& m, int i, int j, int k);
    int hterminal(Molecule& m, BitArrayLTri& bonds, int i);
    int nearest_contact(int i, Molecule& m);

    void add_bonds(const Ref<SetIntCoor>& list, BitArrayLTri& bonds, Molecule& m);
    void add_bends(const Ref<SetIntCoor>& list, BitArrayLTri& bonds, Molecule& m);
    void add_tors(const Ref<SetIntCoor>& list, BitArrayLTri& bonds, Molecule& m);
    void add_out(const Ref<SetIntCoor>& list, BitArrayLTri& bonds, Molecule& m);

    /// (potentially) modifies the adjacency matrix to make sure that there are no disconnected subgraphs
    void connect_subgraphs(const Molecule& mol,
                           BitArrayLTri& adjacency_matrix);

  public:
    /** Create an IntCoorGen given a Molecule and, optionally, extra bonds.
        IntCoorGen keeps a reference to extra and deletes it when the
        destructor is called. */
    IntCoorGen(const Ref<Molecule>&, int nextra=0, int *extra=0);
    /** The KeyVal constructor.
        <dl>

        <dt><tt>molecule</tt><dd> A Molecule object.  There is no default.

        <dt><tt>radius_scale_factor</tt><dd> If the distance between two
        atoms is less than the radius scale factor times the sum of the
        atoms' atomic radii, then a bond is placed between the two atoms
        for the purpose of finding internal coordinates.  The default is
        1.1.

        <dt><tt>linear_bend_threshold</tt><dd> A bend angle in degress
        greater than 180 minus this keyword's floating point value is
        considered a linear bend. The default is 1.0.

        <dt><tt>linear_tors_threshold</tt><dd> The angles formed by atoms
        a-b-c and b-c-d are checked for near linearity.  If an angle in
        degrees is greater than 180 minus this keyword's floating point
        value, then the torsion is classified as a linear torsion. The
        default is 1.0.

        <dt><tt>linear_bend</tt><dd> Generate BendSimpleCo objects to
        describe linear bends.  The default is false.

        <dt><tt>linear_lbend</tt><dd> Generate pairs of LinIPSimpleCo and
        LinIPSimpleCo objects to describe linear bends.  The default is
        true.

        <dt><tt>linear_tors</tt><dd> Generate TorsSimpleCo objects to
        described linear torsions.  The default is false.

        <dt><tt>linear_stors</tt><dd> Generate ScaledTorsSimpleCo objects
        to described linear torsions.  The default is true.

        <dt><tt>extra_bonds</tt><dd> This is a vector of atom numbers,
        where elements \f$2 (i-1) + 1\f$ and \f$2 i\f$ specify the atoms
        which are bound in extra bond \f$i\f$.  The extra_bonds keyword
        should only be needed for weakly interacting fragments, otherwise
        all the needed bonds will be found.

        </dl> */
    IntCoorGen(const Ref<KeyVal>&);
    IntCoorGen(StateIn&);

    ~IntCoorGen();

    /// Standard member.
    void save_data_state(StateOut&);

    /// This generates a set of internal coordinates.
    virtual void generate(const Ref<SetIntCoor>&);

    /// Print out information about this.
    virtual void print(std::ostream& out=ExEnv::out0()) const;

    /// computes the adjacency matrix for this molecule using atomic radii
    /// and the scaling_factor
    static BitArrayLTri adjacency_matrix(const Molecule& mol,
                                         double radius_scaling_factor = 1.1);
    /// given the adjacency matrix find all disconnected subgraphs,
    /// each subgraph is specified by a set of vertex indices
    static std::vector<std::set<int> > find_disconnected_subgraphs(const BitArrayLTri& adjmat);
};


// ////////////////////////////////////////////////////////////////////////


/** The MolecularCoor abstract class describes the coordinate system used
to describe a molecule.  It is used to convert a molecule's cartesian
coordinates to and from this coordinate system. */
class MolecularCoor: public SavableState
{
  protected:
    Ref<Molecule> molecule_;
    RefSCDimension dnatom3_; // the number of atoms x 3
    Ref<SCMatrixKit> matrixkit_; // used to construct matrices

    int debug_;
  public:
    MolecularCoor(Ref<Molecule>&);
    MolecularCoor(StateIn&);
    /** The KeyVal constructor.
        <dl>

        <dt><tt>molecule</tt><dd> A Molecule object.  There is no default.

        <dt><tt>debug</tt><dd> An integer which, if nonzero, will cause
        extra output.

        <dt><tt>matrixkit</tt><dd> A SCMatrixKit object.  It is usually
        unnecessary to give this keyword.

        <dt><tt>natom3</tt><dd> An SCDimension object for the dimension of
        the vector of cartesian coordinates.  It is usually unnecessary to
        give this keyword.

        </dl> */
    MolecularCoor(const Ref<KeyVal>&);

    virtual ~MolecularCoor();

    void save_data_state(StateOut&);

    /** Returns a smart reference to an SCDimension equal to the
        number of atoms in the molecule times 3. */
    RefSCDimension dim_natom3() { return dnatom3_; }

    /// Returns the molecule.
    Ref<Molecule> molecule() const { return molecule_; }

    /// Print the coordinate.
    virtual void print(std::ostream& =ExEnv::out0()) const = 0;
    virtual void print_simples(std::ostream& =ExEnv::out0()) const = 0;

    /** Returns a smart reference to an SCDimension equal to the number of
        coordinates (be they Cartesian, internal, or whatever) that are
        being optimized. */
    virtual RefSCDimension dim() = 0;
    
    /** Given a set of displaced internal coordinates, update the cartesian
        coordinates of the Molecule contained herein.  This function does
        not change the vector ``internal''. */
    int to_cartesian(const RefSCVector&internal);
    virtual int to_cartesian(const Ref<Molecule>&mol,
                             const RefSCVector&internal) = 0;

    /** Fill in the vector ``internal'' with the current internal
        coordinates.  Note that this member will update the values of the
        variable internal coordinates. */
    virtual int to_internal(RefSCVector&internal) = 0;

    /** Convert the internal coordinate gradients in ``internal'' to
        Cartesian coordinates and copy these Cartesian coordinate gradients
        to ``cartesian''. Only the variable internal coordinate gradients
        are transformed. */
    virtual int to_cartesian(RefSCVector&cartesian,RefSCVector&internal) = 0;

    /** Convert the Cartesian coordinate gradients in ``cartesian'' to
        internal coordinates and copy these internal coordinate gradients
        to ``internal''.  Only the variable internal coordinate gradients
        are calculated. */
    virtual int to_internal(RefSCVector&internal,RefSCVector&cartesian) = 0;

    /** Convert the internal coordinate Hessian ``internal'' to Cartesian
        coordinates and copy the result to ``cartesian''.  Only the variable
        internal coordinate force constants are transformed. */
    virtual int to_cartesian(RefSymmSCMatrix&cartesian,
                              RefSymmSCMatrix&internal) =0;

    /** Convert the Cartesian coordinate Hessian ``cartesian'' to internal
        coordinates and copy the result to ``internal''.  Only the variable
        internal coordinate force constants are calculated. */
    virtual int to_internal(RefSymmSCMatrix&internal,
                             RefSymmSCMatrix&cartesian) = 0;

    /** Calculate an approximate hessian and place the result in
        ``hessian''. */
    virtual void guess_hessian(RefSymmSCMatrix&hessian) = 0;

    /** Given an Hessian, return the inverse of that hessian.  For singular
        matrices this should return the generalized inverse. */
    virtual RefSymmSCMatrix inverse_hessian(RefSymmSCMatrix&) = 0;

    /// Returns the number of constrained coordinates.
    virtual int nconstrained();

    /** When this is called, MoleculeCoor may select a new internal
        coordinate system and return a transform to it.  The default action
        is to not change anything and return an IdentityTransform. */
    virtual Ref<NonlinearTransform> change_coordinates();

    Ref<SCMatrixKit> matrixkit() const { return matrixkit_; }
};


/** The IntMolecularCoor abstract class describes a molecule's coordinates
in terms of internal coordinates. */
class IntMolecularCoor: public MolecularCoor
{
  protected:
    Ref<IntCoorGen> generator_;

    void form_K_matrix(RefSCDimension& dredundant,
                       RefSCDimension& dfixed,
                       RefSCMatrix& K,
                       int*& is_totally_symmetric);

    RefSCDimension dim_; // corresponds to the number of variable coordinates
    RefSCDimension dvc_; // the number of variable + constant coordinates

    Ref<SetIntCoor> variable_; // the variable internal coordinates
    Ref<SetIntCoor> constant_; // the constant internal coordinates
    
    Ref<SetIntCoor> fixed_;
    Ref<SetIntCoor> watched_;
    Ref<IntCoor> followed_;

    // these are all of the basic coordinates
    Ref<SetIntCoor> bonds_;
    Ref<SetIntCoor> bends_;
    Ref<SetIntCoor> tors_;
    Ref<SetIntCoor> outs_;
    // these are provided by the user or generated coordinates that
    // could not be assigned to any of the above catagories
    Ref<SetIntCoor> extras_;

    Ref<SetIntCoor> all_;

    // Useful relationships
    // variable_->n() + constant_->n() = 3N-6(5)
    // symm_->n() + asymm_->n() = 3N-6(5)

    int update_bmat_;  // if 1 recompute the b matrix during to_cartesian
    int only_totally_symmetric_; // only coors with tot. symm comp. are varied
    double symmetry_tolerance_; // tol used to find coors with tot. sym. comp.
    double simple_tolerance_; // tol used to see if a simple is included
    double coordinate_tolerance_; // tol used to see if a coor is included
    double cartesian_tolerance_;  // tol used in intco->cart transformation
    double scale_bonds_; // scale factor for bonds
    double scale_bends_; // scale factor for bends
    double scale_tors_;  // scale factor for tors
    double scale_outs_;  // scale factor for outs

    int nextra_bonds_;
    int* extra_bonds_;

    int given_fixed_values_; // if true make molecule have given fixed values

    int decouple_bonds_;
    int decouple_bends_;

    int max_update_steps_;
    double max_update_disp_;

    /** This is called by the constructors of classes derived from
        IntMolecularCoor.  It initialized the lists of simple internal
        coordinates, and then calls the form_coordinates() member. */
    virtual void init();
    /** Allocates memory for the SetIntCoor's used to store the
        simple and internal coordinates. */
    virtual void new_coords();
    /// Reads the KeyVal input.
    virtual void read_keyval(const Ref<KeyVal>&);

    // control whether or not to print coordinates when they are formed
    int form_print_simples_;
    int form_print_variable_;
    int form_print_constant_;
    int form_print_molecule_;
  public:
    IntMolecularCoor(StateIn&);
    IntMolecularCoor(Ref<Molecule>&mol);
    /** The KeyVal constructor.
        <dl>

        <dt><tt>variable</tt><dd> Gives a SetIntCoor object that specifies
        the internal coordinates that can be varied. If this is not given,
        the variable coordinates will be generated.

        <dt><tt>followed</tt><dd> Gives a IntCoor object that specifies a
        coordinate to used as the first coordinate in the variable
        coordinate list.  The remaining coordinates will be automatically
        generated.  The default is no followed coordinate.  This option is
        usually used to set the initial search direction for a transition
        state optimization, where it is used in conjunction with the
        mode_following keyword read by the EFCOpt class.

        <dt><tt>fixed</tt><dd> Gives a SetIntCoor object that specifies the
        internal coordinates that will be fixed.  The default is no fixed
        coordinates.

        <dt><tt>watched</tt><dd> Gives a SetIntCoor object that specifies
        internal coordinates that will be printed out whenever the
        coordinates are changed.  The default is none.

        <dt><tt>have_fixed_values</tt><dd> If true, then values for the
        fixed coordinates must be given in fixed and an attempt will be
        made to displace the initial geometry to the given fixed
        values. The default is false.

        <dt><tt>extra_bonds</tt><dd> This is only read if the generator
        keyword is not given.  It is a vector of atom numbers, where
        elements \f$(i-1)\times 2 + 1\f$ and \f$i\times 2\f$ specify the
        atoms which are bound in extra bond \f$i\f$.  The extra_bonds
        keyword should only be needed for weakly interacting fragments,
        otherwise all the needed bonds will be found.

        <dt><tt>generator</tt><dd> Specifies an IntCoorGen object that
        creates simple, redundant internal coordinates. If this keyword is
        not given, then a vector giving extra bonds to be added is read
        from extra_bonds and this is used to create an IntCoorGen object.

        <dt><tt>decouple_bonds</tt><dd> Automatically generated internal
        coordinates are linear combinations of possibly any mix of simple
        internal coordinates.  If decouple_bonds is true, an attempt will
        be made to form some of the internal coordinates from just stretch
        simple coordinates.  The default is false.

        <dt><tt>decouple_bends</tt><dd> This is like decouple_bonds except
        it applies to the bend-like coordinates.  The default is false.

        <dt><tt>max_update_disp</tt><dd> The maximum displacement to be
        used in the displacement to fixed internal coordinates values.
        Larger displacements will be broken into several smaller
        displacements and new coordinates will be formed for each of these
        displacments. This is only used when fixed and have_fixed_values
        are given.  The default is 0.5.

        <dt><tt>max_update_steps</tt><dd> The maximum number of steps
        permitted to convert internal coordinate displacements to cartesian
        coordinate displacements.  The default is 100.
               
        <dt><tt>update_bmat</tt><dd> Displacements in internal coordinates
        are converted to a cartesian displacements iteratively.  If there
        are large changes in the cartesian coordinates during conversion,
        then recompute the \f$B\f$ matrix, which is using to do the
        conversion.  The default is false.

        <dt><tt>only_totally_symmetric</tt><dd> If a simple test reveals
        that an internal coordinate is not totally symmetric, then it will
        not be added to the internal coordinate list.  The default is true.
               
        <dt><tt>simple_tolerance</tt><dd> The internal coordinates are
        formed as linear combinations of simple, redundant internal
        coordinates.  Coordinates with coefficients smaller then
        simple_tolerance will be omitted. The default is 1.0e-3.

        <dt><tt>cartesian_tolerance</tt><dd> The tolerance for conversion
        of internal coordinate displacements to cartesian displacements.
        The default is 1.0e-12.
               
        <dt><tt>form:print_simple</tt><dd> Print the simple internal
        coordinates.  The default is false.
               
        <dt><tt>form:print_variable</tt><dd> Print the variable internal
        coordinates.  The default is false.
               
        <dt><tt>form:print_constant</tt><dd> Print the constant internal
        coordinates.  The default is false.
               
        <dt><tt>form:print_molecule</tt><dd> Print the molecule when
        forming coordinates.  The default is false.
               
        <dt><tt>scale_bonds</tt><dd> Obsolete.  The default value is 1.0.

        <dt><tt>scale_bends</tt><dd> Obsolete.  The default value is 1.0.
               
        <dt><tt>scale_tors</tt><dd> Obsolete.  The default value is 1.0.

        <dt><tt>scale_outs</tt><dd> Obsolete.  The default value is 1.0.

        <dt><tt>symmetry_tolerance</tt><dd> Obsolete.  The default is 1.0e-5.

        <dt><tt>coordinate_tolerance</tt><dd> Obsolete.  The default is 1.0e-7.

        </dl> */
    IntMolecularCoor(const Ref<KeyVal>&);

    virtual ~IntMolecularCoor();
    void save_data_state(StateOut&);
  
    /** Actually form the variable and constant internal coordinates from
        the simple internal coordinates. */
    virtual void form_coordinates(int keep_variable=0) =0;
    
    /** Like to_cartesians(), except all internal coordinates are
        considered, not just the variable ones. */
    virtual int all_to_cartesian(const Ref<Molecule> &,RefSCVector&internal);
    /** Like to_internal(), except all internal coordinates are
        considered, not just the variable ones. */
    virtual int all_to_internal(const Ref<Molecule> &,RefSCVector&internal);

    /** These implement the virtual functions inherited from
        MolecularCoor. */
    virtual RefSCDimension dim();
    virtual int to_cartesian(const Ref<Molecule> &,const RefSCVector&internal);
    virtual int to_internal(RefSCVector&internal);
    virtual int to_cartesian(RefSCVector&cartesian,RefSCVector&internal);
    virtual int to_internal(RefSCVector&internal,RefSCVector&cartesian);
    virtual int to_cartesian(RefSymmSCMatrix&cart,RefSymmSCMatrix&internal);
    virtual int to_internal(RefSymmSCMatrix&internal,RefSymmSCMatrix&cart);
    virtual void print(std::ostream& =ExEnv::out0()) const;
    virtual void print_simples(std::ostream& =ExEnv::out0()) const;
    virtual void print_variable(std::ostream& =ExEnv::out0()) const;
    virtual void print_constant(std::ostream& =ExEnv::out0()) const;
    int nconstrained();
};

// ///////////////////////////////////////////////////////////////////////

/** The SymmMolecularCoor class derives from IntMolecularCoor.  It provides
a unique set of totally symmetric internal coordinates.  Giving an
MolecularEnergy object a coor is usually the best way to optimize a
molecular structure.  However, for some classes of molecules
SymmMolecularCoor doesn't work very well.  For example, enediyne can cause
problems.  In these cases, cartesian coordinates (obtained by not giving
the MolecularEnergy object the coor keyword) might be better or you can
manually specify the coordinates that the SymmMolecularCoor object uses
with the variable keyword (see the IntMolecularCoor class description).  */
class SymmMolecularCoor: public IntMolecularCoor
{
  protected:
    // true if coordinates should be changed during optimization
    int change_coordinates_;
    // true if hessian should be transformed too (should always be true)
    int transform_hessian_;
    // max value for the condition number if coordinates can be changed
    double max_kappa2_;

    void init();
  public:
    SymmMolecularCoor(Ref<Molecule>&mol);
    SymmMolecularCoor(StateIn&);
    /** The KeyVal constructor.
        <dl>

        <dt><tt>change_coordinates</tt><dd> If true, the quality of the
        internal coordinates will be checked periodically and if they are
        beginning to become linearly dependent a new set of internal
        coordinates will be computed.  The default is false.

        <dt><tt>max_kappa2</tt><dd> A measure of the quality of the
        internal coordinates.  Values of the 2-norm condition,
        \f$\kappa_2\f$, larger than max_kappa2 are considered linearly
        dependent.  The default is 10.0.

        <dt><tt>transform_hessian</tt><dd> If true, the hessian will be
        transformed every time the internal coordinates are formed.  The
        default is true.

        </dl> */
    SymmMolecularCoor(const Ref<KeyVal>&);

    virtual ~SymmMolecularCoor();
    void save_data_state(StateOut&);

    /** Actually form the variable and constant internal coordinates from
        simple internal coordinates. */
    void form_coordinates(int keep_variable=0);

    /// Form the approximate hessian.
    void guess_hessian(RefSymmSCMatrix&hessian);
    /// Invert the hessian.
    RefSymmSCMatrix inverse_hessian(RefSymmSCMatrix&);

    /** This overrides MoleculeCoor's change_coordinates
        and might transform to a new set of coordinates. */
    Ref<NonlinearTransform> change_coordinates();

    void print(std::ostream& =ExEnv::out0()) const;
};

// ///////////////////////////////////////////////////////////////////////

/** The RedundMolecularCoor class provides a redundant set of simple
internal coordinates. */
class RedundMolecularCoor: public IntMolecularCoor
{

  public:
    RedundMolecularCoor(Ref<Molecule>&mol);
    RedundMolecularCoor(StateIn&);
    /// The KeyVal constructor.
    RedundMolecularCoor(const Ref<KeyVal>&);

    virtual ~RedundMolecularCoor();
    void save_data_state(StateOut&);

    /** Actually form the variable and constant internal coordinates from
        the simple internal coordinates. */
    void form_coordinates(int keep_variable=0);
    /// Form the approximate hessian.
    void guess_hessian(RefSymmSCMatrix&hessian);
    /// Invert the hessian.
    RefSymmSCMatrix inverse_hessian(RefSymmSCMatrix&);
};

// ///////////////////////////////////////////////////////////////////////

/** The CartMolecularCoor class implements Cartesian coordinates in a way
    suitable for use in geometry optimizations. CartMolecularCoor is a
    SavableState has StateIn and KeyVal constructors. CartMolecularCoor is
    derived from MolecularCoor. */
class CartMolecularCoor: public MolecularCoor
{
  private:
  protected:
    RefSCDimension dim_; // the number of atoms x 3

    /// Initializes the dimensions.
    virtual void init();
  public:
    CartMolecularCoor(Ref<Molecule>&mol);
    CartMolecularCoor(StateIn&);
    /// The KeyVal constructor.
    CartMolecularCoor(const Ref<KeyVal>&);

    virtual ~CartMolecularCoor();
  
    void save_data_state(StateOut&);

    /// These implement the virtual functions inherited from MolecularCoor.
    virtual RefSCDimension dim();
    virtual int to_cartesian(const Ref<Molecule>&,const RefSCVector&internal);
    virtual int to_internal(RefSCVector&internal);
    virtual int to_cartesian(RefSCVector&cartesian,RefSCVector&internal);
    virtual int to_internal(RefSCVector&internal,RefSCVector&cartesian);
    virtual int to_cartesian(RefSymmSCMatrix&cart,RefSymmSCMatrix&internal);
    virtual int to_internal(RefSymmSCMatrix&internal,RefSymmSCMatrix&cart);
    virtual void print(std::ostream& =ExEnv::out0()) const;
    virtual void print_simples(std::ostream& =ExEnv::out0()) const;
    void guess_hessian(RefSymmSCMatrix&hessian);
    RefSymmSCMatrix inverse_hessian(RefSymmSCMatrix&);
};

/// @}
// end of addtogroup ChemistryMolecule

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
