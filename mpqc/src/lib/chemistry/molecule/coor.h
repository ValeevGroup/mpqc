
#ifndef _chemistry_molecule_coor_h
#define _chemistry_molecule_coor_h

#ifdef __GNUC__
#pragma interface
#include <ostream.h>
#endif

#include <util/container/array.h>
#include <util/container/set.h>
#include <math/scmat/matrix.h>
#include <math/optimize/transform.h>
#include <chemistry/molecule/molecule.h>

class SSRefIntCoor;
typedef class SSRefIntCoor RefIntCoor;

//.  \clsnm{IntCoor} is an abstract base class.  From it are derived the
//simple internal coordinate classes (derivatives of \clsnmref{SimpleCo})
//and a class describing linear combinations of internal coordinates
//(\clsnmref{SumIntCoor}).
//
// \clsnm{IntCoor} is a \clsnmref{SavableState} and has a
//clsnmref{StateIn} constructor, as well as a \clsnmref{KeyVal}
//constructor.
class IntCoor: public SavableState {
#   define CLASSNAME IntCoor
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    // conversion factors from radians, bohr to the preferred units
    static double bohr_conv;
    static double radian_conv;
    char *label_;
    double value_;
  public:
    IntCoor(StateIn&);
    IntCoor(const IntCoor&);
    //. This constructor takes a string containing a label for the
    // internal coordinate.  The string is copied.
    IntCoor(const char* label = 0);
    //. The \clsnmref{KeyVal} constructor.
    IntCoor(const RefKeyVal&);
    
    virtual ~IntCoor();
    void save_data_state(StateOut&);

    //. Returns the string containing the label for the internal coordinate.
    virtual const char* label() const;
    //. Returns the value of the coordinate in atomic units or radians.
    virtual double value() const;
    //. Returns the value of the coordinate in more familiar units (at
    // least to those in the U.S.).
    virtual double preferred_value() const;
    //. Returns a string representation of the type of coordinate this is.
    virtual const char* ctype() const = 0;
#ifdef __GNUC__
    //. Print information about the coordinate.
    virtual void print(RefMolecule =0, ostream& =cout);
#else
    virtual void print();
    virtual void print(RefMolecule, ostream& =cout);
#endif
    //. Returns the value of the force constant associated with this
    // coordinate.
    virtual double force_constant(RefMolecule&) = 0;
    //. Recalculate the value of the coordinate.
    virtual void update_value(const RefMolecule&) = 0;
    //. Fill in a row the the B matrix.
    virtual void bmat(RefMolecule&,RefSCVector&bmat,double coef = 1.0) = 0;
    //. Test to see if this internal coordinate is equivalent to that one.
    // The definition of equivalence is left up to the individual coordinates.
    virtual int equivalent(RefIntCoor&) = 0;
};
SavableState_REF_dec(IntCoor);
ARRAY_dec(RefIntCoor);
SET_dec(RefIntCoor);
ARRAYSET_dec(RefIntCoor);

//.  \clsnm{SumIntCoor} is used to construct linear combinations of
//internal coordinates.  Normally one will use simple internal coordinates,
//such as bond lengths and angles.  \clsnm{SumIntCoor} inherits from
//\clsnmref{IntCoor}, so it is a \clsnmref{SavableState}.
//\clsnm{SumIntCoor} has \clsnmref{StateIn} and \clsnmref{KeyVal}
//constructors.
class SumIntCoor: public IntCoor {
#   define CLASSNAME SumIntCoor
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    Arraydouble coef_;
    ArrayRefIntCoor coor_;
  public:
    SumIntCoor(StateIn&);
    //. This constructor takes a string containing a label for this
    // coordinate.
    SumIntCoor(const char *);
    //. The \clsnmref{KeyVal} constructor.
    SumIntCoor(const RefKeyVal&);

    ~SumIntCoor();
    void save_data_state(StateOut&);

    //. Returns the number of coordinates in this linear combination.
    int n();
    //. Add a coordinate to the linear combination.  \vrbl{coef} is the
    // coefficient for the added coordinate.
    void add(RefIntCoor&,double coef);
    //. This function normalizes all the coefficients.
    void normalize();

    // IntCoor overrides
    //. Returns the value of the coordinate in a.u. and radians.
    double preferred_value() const;
    //. Always returns ``SUM''.
    const char* ctype() const;
#ifdef __GNUC__
    //. Print the individual coordinates in the sum with their coefficients.
    void print(RefMolecule = 0, ostream& =cout);
#else
    void print();
    void print(RefMolecule, ostream& =cout);
#endif
    //. Returns the weighted sum of the individual force constants.
    double force_constant(RefMolecule&);
    //. Recalculate the value of the coordinate.
    void update_value(const RefMolecule&);
    //. Fill in a row the the B matrix.
    void bmat(RefMolecule&,RefSCVector&bmat,double coef = 1.0);
    //. Always returns 0.
    int equivalent(RefIntCoor&);
};

class SSRefSetIntCoor;
typedef class SSRefSetIntCoor RefSetIntCoor;

//.  \clsnm{SetIntCoor} is a class which holds sets of internal
//coordinates, be they simple internal coordinates or combinations of
//coordinates.  \clsnm{SetIntCoor} is a \clsnmref{SavableState}, and has
//\clsnmref{StateIn} and \clsnmref{KeyVal} constructors.
class SetIntCoor: public SavableState {
#   define CLASSNAME SetIntCoor
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    ArraysetRefIntCoor coor_;
  public:
    SetIntCoor();
    SetIntCoor(StateIn&);
    //. The \clsnmref{KeyVal} constructor.
    SetIntCoor(const RefKeyVal&);

    virtual ~SetIntCoor();
    void save_data_state(StateOut&);

    //. Adds an internal coordinate to the set.
    void add(const RefIntCoor&);
    //. Adds all the elements of another set to this one.
    void add(const RefSetIntCoor&);
    //. Removes a coordinate from this set.
    void del(const RefIntCoor&);
    //. Removes all the elements of a set of coordinates from this one.
    void del(const RefSetIntCoor&);
    //. Removes all coordinates from the set.
    void clear();
    //. Returns the number of coordinates in the set.
    int n() const;
    //. Returns a reference to the i'th coordinate in the set.
    RefIntCoor coor(int i) const;
    //. Compute the B matrix by finite displacements.
    virtual void fd_bmat(RefMolecule&,RefSCMatrix&);
    //. Compute the B matrix the old-fashioned way.
    virtual void bmat(RefMolecule&, RefSCMatrix&);
    //. Create an approximate Hessian for this set of coordinates.  This
    // Hessian is a symmetric matrix whose i'th diagonal is the force constant
    // for the i'th coordinate in the set.
    virtual void guess_hessian(RefMolecule&,RefSymmSCMatrix&);
#ifdef __GNUC__
    //. Print the coordinates in the set.
    virtual void print(RefMolecule =0,ostream& =cout);
#else
    virtual void print();
    virtual void print(RefMolecule,ostream& =cout);
#endif
    //. Recalculate the values of the internal coordinates in the set.
    virtual void update_values(const RefMolecule&);
    //. Copy the values of the internal coordinates to a vector.
    virtual void values_to_vector(const RefSCVector&);
};
SavableState_REF_dec(SetIntCoor);

//////////////////////////////////////////////////////////////////////////

class BitArray;

//.  \clsnm{IntCoorGen} generates a set of simple internal coordinates
//given a molecule.
class IntCoorGen: public SavableState
{
#   define CLASSNAME IntCoorGen
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    RefMolecule molecule_;
    RefMessageGrp grp_;
    
    int linear_bends_;
    int linear_tors_;
    int linear_stors_;
    int nextra_bonds_;
    int *extra_bonds_;
    double linear_bend_thres_;
    double linear_tors_thres_;
    double radius_scale_factor_;

    double cos_ijk(Molecule& m, int i, int j, int k);
    int hterminal(Molecule& m, BitArray& bonds, int i);
    int nearest_contact(int i, Molecule& m);

    void add_bonds(const RefSetIntCoor& list, BitArray& bonds, Molecule& m);
    void add_bends(const RefSetIntCoor& list, BitArray& bonds, Molecule& m);
    void add_tors(const RefSetIntCoor& list, BitArray& bonds, Molecule& m);
    void add_out(const RefSetIntCoor& list, BitArray& bonds, Molecule& m);
  public:
    //. Create an \clsnm{IntCoorGen} given a \clsnmref{Molecule} and,
    //optionally, extra bonds.  \clsnm{IntCoorGen} keeps a reference to
    //\vrbl{extra} and deletes it when the destructor is called.
    IntCoorGen(const RefMolecule&, int nextra=0, int *extra=0);
    //. Standard constructors for \clsnm{IntCoorGen}.
    IntCoorGen(const RefKeyVal&);
    IntCoorGen(StateIn&);

    ~IntCoorGen();

    //. Standard member.
    void save_data_state(StateOut&);

    //. This generates a set of internal coordinates.
    virtual void generate(const RefSetIntCoor&);

    //. Set the message group used by the coordinate generator
    void set_messagegrp(const RefMessageGrp& g) { grp_=g; }
    
    //. Print out information about this.
    virtual void print(ostream& out=cout);
};
SavableState_REF_dec(IntCoorGen);

//////////////////////////////////////////////////////////////////////////

//.  \clsnm{MolecularCoor} is a virtual base class for a set of classes
//used to describe coordinates for molecules useful in optimizing
//geometries (among other things).  \clsnm{MolecularCoor} is a
//\clsnmref{SavableState} and has \clsnmref{StateIn} and \clsnmref{KeyVal}
//constructors.
class MolecularCoor: public SavableState
{
#   define CLASSNAME MolecularCoor
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    RefMolecule molecule_;
    RefSCDimension dnatom3_; // the number of atoms x 3
    RefSCMatrixKit matrixkit_; // used to construct matrices
  public:
    MolecularCoor(RefMolecule&);
    MolecularCoor(StateIn&);

    //. The \clsnmref{KeyVal} constructor.
    MolecularCoor(const RefKeyVal&);

    virtual ~MolecularCoor();

    void save_data_state(StateOut&);

    //. Returns a smart reference to an \clsnmref{SCDimension} equal to the
    // number of atoms in the molecule times 3.
    RefSCDimension dim_natom3() { return dnatom3_; }

    //. Print the coordinate.
    virtual void print(ostream& =cout) = 0;
    virtual void print_simples(ostream& =cout) = 0;

    //. Returns a smart reference to an \clsnmref{SCDimension} equal to the
    // number of coordinates (be they Cartesian, internal, or whatever)
    // that are being optimized.
    virtual RefSCDimension dim() = 0;
    
    //. Given a set of displaced internal coordinates, update the cartesian
    // coordinates of the \clsnmref{Molecule} contained herein.  This function
    // does not change the vector ``internal''.
    virtual int to_cartesian(const RefSCVector&internal) = 0;

    //. Fill in the vector ``internal'' with the current internal
    // coordinates.  Note that this member will update the values of the
    // variable internal coordinates.
    virtual int to_internal(RefSCVector&internal) = 0;

    //. Convert the Cartesian coordinate gradients in ``cartesian'' to
    // internal coordinates and copy these internal coordinate gradients
    // to ``internal''.  Only the variable internal coordinate gradients
    // are calculated.
    virtual int to_cartesian(RefSCVector&cartesian,RefSCVector&internal) = 0;

    //. Convert the internal coordinate gradients in ``internal'' to
    // Cartesian coordinates and copy these Cartesian coordinate gradients
    // to ``cartesian''. Only the variable internal coordinate gradients
    // are transformed.
    virtual int to_internal(RefSCVector&internal,RefSCVector&cartesian) = 0;

    //. Convert the internal coordinate Hessian ``internal'' to Cartesian
    // coordinates and copy the result to ``cartesian''.  Only the variable
    // internal coordinate force constants are transformed.
    virtual int to_cartesian(RefSymmSCMatrix&cartesian,
                              RefSymmSCMatrix&internal) =0;

    //. Convert the Cartesian coordinate Hessian ``cartesian'' to internal
    // coordinates and copy the result to ``internal''.  Only the variable
    // internal coordinate force constants are calculated.
    virtual int to_internal(RefSymmSCMatrix&internal,
                             RefSymmSCMatrix&cartesian) = 0;

    //. Calculate an approximate hessian and place the result in
    // ``hessian''.
    virtual void guess_hessian(RefSymmSCMatrix&hessian) = 0;

    //. Given an Hessian, return the inverse of that hessian.  For singular
    // matrices this should return the generalized inverse.
    virtual RefSymmSCMatrix inverse_hessian(RefSymmSCMatrix&) = 0;

    //. Returns the number of constrained coordinates.
    virtual int nconstrained();

    //. When this is called, \clsnmref{MoleculeCoor} may select a new
    // internal coordinate system and return a transform to it.
    // The default action is
    // to not change anything and return an \clsnmref{IdentityTransform}.
    virtual RefNonlinearTransform change_coordinates();

    RefSCMatrixKit matrixkit() { return matrixkit_; }
};
SavableState_REF_dec(MolecularCoor);

//.  \clsnm{IntMolecularCoor} is a virtual base class for the internal
//coordinate classes.  Internal coordinates are very useful in geometry
//optimizations, and are also used to describe the normal vibrational modes
//in spectroscopy.  \clsnm{IntMolecularCoor} is a \clsnmref{SavableState}
//and has \clsnmref{StateIn} and \clsnmref{KeyVal} constructors.
class IntMolecularCoor: public MolecularCoor
{
#   define CLASSNAME IntMolecularCoor
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    RefIntCoorGen generator_;

    void form_K_matrices(RefSCDimension& dredundant,
                         RefSCDimension& dfixed,
                         RefSCMatrix& K,
                         RefSCMatrix& Kfixed,
                         int*& is_totally_symmetric);

    RefSCDimension dim_; // corresponds to the number of variable coordinates
    RefSCDimension dvc_; // the number of variable + constant coordinates

    RefSetIntCoor variable_; // the variable internal coordinates
    RefSetIntCoor constant_; // the constant internal coordinates
    
    RefSetIntCoor fixed_;
    RefIntCoor followed_;

    // these are all of the basic coordinates
    RefSetIntCoor bonds_;
    RefSetIntCoor bends_;
    RefSetIntCoor tors_;
    RefSetIntCoor outs_;
    // these are provided by the user or generated coordinates that
    // could not be assigned to any of the above catagories
    RefSetIntCoor extras_;

    RefSetIntCoor all_;

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

    //. This is called by the constructors of classes derived from
    // \clsnm{IntMolecularCoor}.  It initialized the lists of simple internal
    // coordinates, and then calls the \srccd{form\_coordinates()} member.
    virtual void init();
    //. Allocates memory for the \clsnmref{SetIntCoor}s used to store the
    // simple and internal coordinates.
    virtual void new_coords();
    //. Reads the \clsnmref{KeyVal} input.
    virtual void read_keyval(const RefKeyVal&);
  public:
    IntMolecularCoor(StateIn&);
    IntMolecularCoor(RefMolecule&mol);
    //. The \clsnmref{KeyVal} constructor.
    IntMolecularCoor(const RefKeyVal&);

    virtual ~IntMolecularCoor();
    void save_data_state(StateOut&);
  
    //. Actually form the variable and constant internal coordinates from
    // the simple internal coordinates.
    virtual void form_coordinates() =0;
    
    //. Like \srccd{to\_cartesians()}, except all internal coordinates are
    // considered, not just the variable ones.
    virtual int all_to_cartesian(RefSCVector&internal);
    //. Like \srccd{to\_internal()}, except all internal coordinates are
    // considered, not just the variable ones.
    virtual int all_to_internal(RefSCVector&internal);

    //. These implement the virtual functions inherited from
    // \clsnmref{MolecularCoor}.
    virtual RefSCDimension dim();
    virtual int to_cartesian(const RefSCVector&internal);
    virtual int to_internal(RefSCVector&internal);
    virtual int to_cartesian(RefSCVector&cartesian,RefSCVector&internal);
    virtual int to_internal(RefSCVector&internal,RefSCVector&cartesian);
    virtual int to_cartesian(RefSymmSCMatrix&cart,RefSymmSCMatrix&internal);
    virtual int to_internal(RefSymmSCMatrix&internal,RefSymmSCMatrix&cart);
    virtual void print(ostream& =cout);
    virtual void print_simples(ostream& =cout);
    int nconstrained();
};

/////////////////////////////////////////////////////////////////////////

//.  The \clsnm{SymmMolecularCoor} class implements symmetry adapted linear
//combinations of internal coordinates. \clsnm{SymmMolecularCoor} is a
//\clsnmref{SavableState} and has \clsnmref{StateIn} and \clsnmref{KeyVal}
//constructors.  \clsnm{SymmMolecularCoor} is derived from
//\clsnmref{IntMolecularCoor}.
class SymmMolecularCoor: public IntMolecularCoor
{
#   define CLASSNAME SymmMolecularCoor
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    // true if coordinates should be changed during optimization
    int change_coordinates_;
    // true if hessian should be transformed too (should always be true)
    int transform_hessian_;
    // max value for the condition number if coordinates can be changed
    double max_kappa2_;

    void init();
  public:
    SymmMolecularCoor(RefMolecule&mol);
    SymmMolecularCoor(StateIn&);
    //. The \clsnmref{KeyVal} constructor.
    SymmMolecularCoor(const RefKeyVal&);

    virtual ~SymmMolecularCoor();
    void save_data_state(StateOut&);

    //. Actually form the variable and constant internal coordinates from
    // the simple internal coordinates.
    void form_coordinates();

    //. Form the approximate hessian.
    void guess_hessian(RefSymmSCMatrix&hessian);
    //. Invert the hessian.
    RefSymmSCMatrix inverse_hessian(RefSymmSCMatrix&);

    //. This overrides \clsnmref{MoleculeCoor}'s \srccd{change\_coordinates}
    // and might transform to a new set of coordinates.
    RefNonlinearTransform change_coordinates();

    void print(ostream& =cout);
};

/////////////////////////////////////////////////////////////////////////

//.  The \clsnm{RedundMolecularCoor} class implements redundant sets of
//internal coordinates. \clsnm{RedundMolecularCoor} is a
//\clsnmref{SavableState} and has \clsnmref{StateIn} and \clsnmref{KeyVal}
//constructors.  \clsnm{RedundMolecularCoor} is derived from
//\clsnmref{IntMolecularCoor}.
class RedundMolecularCoor: public IntMolecularCoor
{
#   define CLASSNAME RedundMolecularCoor
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>

  public:
    RedundMolecularCoor(RefMolecule&mol);
    RedundMolecularCoor(StateIn&);
    //. The \clsnmref{KeyVal} constructor.
    RedundMolecularCoor(const RefKeyVal&);

    virtual ~RedundMolecularCoor();
    void save_data_state(StateOut&);

    //. Actually form the variable and constant internal coordinates from
    // the simple internal coordinates.
    void form_coordinates();
    //. Form the approximate hessian.
    void guess_hessian(RefSymmSCMatrix&hessian);
    //. Invert the hessian.
    RefSymmSCMatrix inverse_hessian(RefSymmSCMatrix&);
};

/////////////////////////////////////////////////////////////////////////

//.  The \clsnm{CartMolecularCoor} class implements Cartesian coordinates
//in a way suitable for use in geometry
//optimizations. \clsnm{CartMolecularCoor} is a \clsnmref{SavableState} and
//has \clsnmref{StateIn} and \clsnmref{KeyVal}
//constructors. \clsnm{CartMolecularCoor} is derived from
//\clsnmref{MolecularCoor}.
class CartMolecularCoor: public MolecularCoor
{
#   define CLASSNAME CartMolecularCoor
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
  protected:
    RefSCDimension dim_; // the number of atoms x 3

    //. Initializes the dimensions.
    virtual void init();
  public:
    CartMolecularCoor(RefMolecule&mol);
    CartMolecularCoor(StateIn&);
    //. The \clsnmref{KeyVal} constructor.
    CartMolecularCoor(const RefKeyVal&);

    virtual ~CartMolecularCoor();
  
    void save_data_state(StateOut&);

    //. These implement the virtual functions inherited from
    // \clsnmref{MolecularCoor}.
    virtual RefSCDimension dim();
    virtual int to_cartesian(const RefSCVector&internal);
    virtual int to_internal(RefSCVector&internal);
    virtual int to_cartesian(RefSCVector&cartesian,RefSCVector&internal);
    virtual int to_internal(RefSCVector&internal,RefSCVector&cartesian);
    virtual int to_cartesian(RefSymmSCMatrix&cart,RefSymmSCMatrix&internal);
    virtual int to_internal(RefSymmSCMatrix&internal,RefSymmSCMatrix&cart);
    virtual void print(ostream& =cout);
    virtual void print_simples(ostream& =cout);
    void guess_hessian(RefSymmSCMatrix&hessian);
    RefSymmSCMatrix inverse_hessian(RefSymmSCMatrix&);
};

#endif
