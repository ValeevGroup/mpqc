
#ifndef _chemistry_molecule_coor_h
#define _chemistry_molecule_coor_h

#ifdef __GNUC__
#pragma interface
#include <ostream.h>
#endif

#include <util/misc/scostream.h>
#include <util/container/ref.h>
#include <util/container/array.h>
#include <util/container/set.h>
#include <math/scmat/matrix.h>
#include <chemistry/molecule/molecule.h>

class SSRefIntCoor;
typedef class SSRefIntCoor RefIntCoor;

//texi
// @code{IntCoor} is an abstract base class.  From it are derived the
// simple internal coordinate classes (@ref{The SimpleCo Class} and
// @ref{The Simple Internal Coordinate Classes}), and a class describing
// linear combinations of internal coordinates (@ref{The SumIntCoor Class}).
//
// @code{IntCoor} is a @code{SavableState} and has a @code{StateIn}
// constructor, as well as a @code{KeyVal} constructor.
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
    //texi This constructor takes a string containing a label for the
    // internal coordinate.  The string is copied.
    IntCoor(const char* label = 0);
    //texi The KeyVal constructor (@ref{The IntCoor KeyVal Constructor}).
    IntCoor(const RefKeyVal&);
    
    virtual ~IntCoor();
    void save_data_state(StateOut&);

    //texi Returns the string containing the label for the internal coordinate.
    virtual const char* label() const;
    //texi Returns the value of the coordinate in atomic units or radians.
    virtual double value() const;
    //texi Returns the value of the coordinate in more familiar units (at
    // least to those in the U.S.).
    virtual double preferred_value() const;
    //texi Returns a string representation of the type of coordinate this is.
    virtual const char* ctype() const = 0;
#ifdef __GNUC__
    //texi Print information about the coordinate.
    virtual void print(RefMolecule =0, SCostream& =SCostream::cout);
#else
    virtual void print();
    virtual void print(RefMolecule, SCostream& =SCostream::cout);
#endif
    //texi Returns the value of the force constant associated with this
    // coordinate.
    virtual double force_constant(RefMolecule&) = 0;
    //texi Recalculate the value of the coordinate.
    virtual void update_value(RefMolecule&) = 0;
    //texi Fill in a row the the B matrix.
    virtual void bmat(RefMolecule&,RefSCVector&bmat,double coef = 1.0) = 0;
    //texi Test to see if this internal coordinate is equivalent to that one.
    // The definition of equivalence is left up to the individual coordinates.
    virtual int equivalent(RefIntCoor&) = 0;
};
SavableState_REF_dec(IntCoor);
ARRAY_dec(RefIntCoor);
SET_dec(RefIntCoor);
ARRAYSET_dec(RefIntCoor);

//texi
// @code{SumIntCoor} is used to construct linear combinations of internal
// coordinates.  Normally one will use simple internal coordinates, such as
// bond lengths and angles.  @code{SumIntCoor} inherits from @code{IntCoor}
// (@ref{The IntCoor Class}), so it is a @code{SavableState}.
// @code{SumIntCoor} has @code{StateIn} and @code{KeyVal} constructors.
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
    //texi This constructor takes a string containing a label for this
    // coordinate.
    SumIntCoor(const char *);
    //texi The KeyVal constructor (@ref{The SumIntCoor KeyVal Constructor}).
    SumIntCoor(const RefKeyVal&);

    ~SumIntCoor();
    void save_data_state(StateOut&);

    //texi Returns the number of coordinates in this linear combination.
    int n();
    //texi Add a coordinate to the linear combination.  @code{coef} is the
    // coefficient for the added coordinate.
    void add(RefIntCoor&,double coef);
    //texi This function normalizes all the coefficients.
    void normalize();

    // IntCoor overrides
    //texi Returns the value of the coordinate in a.u. and radians.
    double preferred_value() const;
    //texi Always returns ``SUM''.
    const char* ctype() const;
#ifdef __GNUC__
    //texi Print the individual coordinates in the sum with their coefficients.
    void print(RefMolecule = 0, SCostream& =SCostream::cout);
#else
    void print();
    void print(RefMolecule, SCostream& =SCostream::cout);
#endif
    //texi Returns the weighted sum of the individual force constants.
    double force_constant(RefMolecule&);
    //texi Recalculate the value of the coordinate.
    void update_value(RefMolecule&);
    //texi Fill in a row the the B matrix.
    void bmat(RefMolecule&,RefSCVector&bmat,double coef = 1.0);
    //texi Always returns 0.
    int equivalent(RefIntCoor&);
};

class SSRefSetIntCoor;
typedef class SSRefSetIntCoor RefSetIntCoor;

//texi
// @code{SetIntCoor} is a class which holds sets of internal coordinates, be
// they simple internal coordinates or combinations of coordinates.
// @code{SetIntCoor} is a @code{SavableState}, and has @code{StateIn} and
// @code{KeyVal} constructors.
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
    //texi The KeyVal constructor (@ref{The SetIntCoor KeyVal Constructor}).
    SetIntCoor(const RefKeyVal&);

    virtual ~SetIntCoor();
    void save_data_state(StateOut&);

    //texi Adds an internal coordinate to the set.
    void add(const RefIntCoor&);
    //texi Adds all the elements of another set to this one.
    void add(const RefSetIntCoor&);
    //texi Removes a coordinate from this set.
    void del(const RefIntCoor&);
    //texi Removes all the elements of a set of coordinates from this one.
    void del(const RefSetIntCoor&);
    //texi Removes all coordinates from the set.
    void clear();
    //texi Returns the number of coordinates in the set.
    int n() const;
    //texi Returns a reference to the i'th coordinate in the set.
    RefIntCoor coor(int i) const;
    //texi Compute the B matrix by finite displacements.
    virtual void fd_bmat(RefMolecule&,RefSCMatrix&);
    //texi Compute the B matrix the old-fashioned way.
    virtual void bmat(RefMolecule&, RefSCMatrix&);
    //texi Create an approximate Hessian for this set of coordinates.  This
    // Hessian is a symmetric matrix whose i'th diagonal is the force constant
    // for the i'th coordinate in the set.
    virtual void guess_hessian(RefMolecule&,RefSymmSCMatrix&);
#ifdef __GNUC__
    //texi Print the coordinates in the set.
    virtual void print(RefMolecule =0,SCostream& =SCostream::cout);
#else
    virtual void print();
    virtual void print(RefMolecule,SCostream& =SCostream::cout);
#endif
    //texi Recalculate the values of the internal coordinates in the set.
    virtual void update_values(RefMolecule&);
    //texi Copy the values of the internal coordinates to a vector.
    virtual void values_to_vector(RefSCVector&);
};
SavableState_REF_dec(SetIntCoor);

//////////////////////////////////////////////////////////////////////////

//texi
// @code{MolecularCoor} is a virtual base class for a set of classes used
// to describe coordinates for molecules useful in optimizing geometries
// (among other things).  @code{MolecularCoor} is a @code{SavableState} and
// has @code{StateIn} and @code{KeyVal} constructors.
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

    //texi The KeyVal constructor (@ref{The MolecularCoor KeyVal Constructor}).
    MolecularCoor(const RefKeyVal&);

    virtual ~MolecularCoor();

    void save_data_state(StateOut&);

    //texi Returns a smart reference to an @code{SCDimension} equal to the
    // number of atoms in the molecule times 3.
    RefSCDimension dim_natom3() { return dnatom3_; }

    //texi Print the coordinate.
    virtual void print(SCostream& =SCostream::cout) = 0;
    virtual void print_simples(SCostream& =SCostream::cout) = 0;

    //texi Returns a smart reference to an @code{SCDimension} equal to the
    // number of coordinates (be they Cartesian, internal, or whatever)
    // that are being optimized.
    virtual RefSCDimension dim() = 0;
    
    //texi Given a set of displaced internal coordinates, update the cartesian
    // coordinates of the @code{Molecule} contained herein.  This function
    // does not change the vector ``internal''.
    virtual int to_cartesian(RefSCVector&internal) = 0;

    //texi Fill in the vector ``internal'' with the current internal
    // coordinates.  Note that this member will update the values of the
    // variable internal coordinates.
    virtual int to_internal(RefSCVector&internal) = 0;

    //texi Convert the Cartesian coordinate gradients in ``cartesian'' to
    // internal coordinates and copy these internal coordinate gradients
    // to ``internal''.  Only the variable internal coordinate gradients
    // are calculated.
    virtual int to_cartesian(RefSCVector&cartesian,RefSCVector&internal) = 0;

    //texi Convert the internal coordinate gradients in ``internal'' to
    // Cartesian coordinates and copy these Cartesian coordinate gradients
    // to ``cartesian''. Only the variable internal coordinate gradients
    // are transformed.
    virtual int to_internal(RefSCVector&internal,RefSCVector&cartesian) = 0;

    //texi Convert the internal coordinate Hessian ``internal'' to Cartesian
    // coordinates and copy the result to ``cartesian''.  Only the variable
    // internal coordinate force constants are transformed.
    virtual int to_cartesian(RefSymmSCMatrix&cartesian,
                              RefSymmSCMatrix&internal) =0;

    //texi Convert the Cartesian coordinate Hessian ``cartesian'' to internal
    // coordinates and copy the result to ``internal''.  Only the variable
    // internal coordinate force constants are calculated.
    virtual int to_internal(RefSymmSCMatrix&internal,
                             RefSymmSCMatrix&cartesian) = 0;

    //texi Calculate an approximate hessian and place the result in
    // ``hessian''.
    virtual void guess_hessian(RefSymmSCMatrix&hessian) = 0;

    //texi Given an Hessian, return the inverse of that hessian.  For singular
    // matrices this should return the generalized inverse.
    virtual RefSymmSCMatrix inverse_hessian(RefSymmSCMatrix&) = 0;
};
SavableState_REF_dec(MolecularCoor);

//texi
// @code{IntMolecularCoor} is a virtual base class for the internal coordinate
// classes.  Internal coordinates are very useful in geometry optimizations,
// and are also used to describe the normal vibrational modes in spectroscopy.
// @code{IntMolecularCoor} is a @code{SavableState} and has @code{StateIn}
// and @code{KeyVal} constructors.
class IntMolecularCoor: public MolecularCoor
{
#   define CLASSNAME IntMolecularCoor
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
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
    RefSetIntCoor extras_; // these are provided by the user

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

    int max_update_steps_;
    double max_update_disp_;

    //texi This is called by the constructors of classes derived from
    // @code{IntMolecularCoor}.  It initialized the lists of simple internal
    // coordinates, and then calls the @code{form_coordinates()} member.
    virtual void init();
    //texi Allocates memory for the @code{SetIntCoor}s used to store the
    // simple and internal coordinates.
    virtual void new_coords();
    //texi Reads the @code{KeyVal} input.
    virtual void read_keyval(const RefKeyVal&);
  public:
    IntMolecularCoor(StateIn&);
    IntMolecularCoor(RefMolecule&mol);
    //texi The KeyVal constructor (@ref{The IntMolecularCoor KeyVal
    // Constructor}).
    IntMolecularCoor(const RefKeyVal&);

    virtual ~IntMolecularCoor();
    void save_data_state(StateOut&);
  
    //texi Actually form the variable and constant internal coordinates from
    // the simple internal coordinates.
    virtual void form_coordinates() =0;
    
    //texi Like @code{to_cartesians()}, except all internal coordinates are
    // considered, not just the variable ones.
    virtual int all_to_cartesian(RefSCVector&internal);
    //texi Like @code{to_internal()}, except all internal coordinates are
    // considered, not just the variable ones.
    virtual int all_to_internal(RefSCVector&internal);

    //texi These implement the virtual functions inherited from
    // @code{MolecularCoor} (@ref{The MolecularCoor Public Interface}).
    virtual RefSCDimension dim();
    virtual int to_cartesian(RefSCVector&internal);
    virtual int to_internal(RefSCVector&internal);
    virtual int to_cartesian(RefSCVector&cartesian,RefSCVector&internal);
    virtual int to_internal(RefSCVector&internal,RefSCVector&cartesian);
    virtual int to_cartesian(RefSymmSCMatrix&cart,RefSymmSCMatrix&internal);
    virtual int to_internal(RefSymmSCMatrix&internal,RefSymmSCMatrix&cart);
    virtual void print(SCostream& =SCostream::cout);
    virtual void print_simples(SCostream& =SCostream::cout);
};

/////////////////////////////////////////////////////////////////////////

//texi
// The @code{SymmMolecularCoor} class implements symmetry adapted linear
// combinations of internal coordinates. @code{SymmMolecularCoor} is a
// @code{SavableState} and has @code{StateIn} and @code{KeyVal} constructors.
// @code{SymmMolecularCoor} is derived from @code{IntMolecularCoor}.
class SymmMolecularCoor: public IntMolecularCoor
{
#   define CLASSNAME SymmMolecularCoor
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>

  public:
    SymmMolecularCoor(RefMolecule&mol);
    SymmMolecularCoor(StateIn&);
    //texi The KeyVal constructor (@ref{The SymmMolecularCoor KeyVal
    // Constructor}).
    SymmMolecularCoor(const RefKeyVal&);

    virtual ~SymmMolecularCoor();
    void save_data_state(StateOut&);

    //texi Actually form the variable and constant internal coordinates from
    // the simple internal coordinates.
    void form_coordinates();

    //texi Form the approximate hessian.
    void guess_hessian(RefSymmSCMatrix&hessian);
    //texi Invert the hessian.
    RefSymmSCMatrix inverse_hessian(RefSymmSCMatrix&);
};

/////////////////////////////////////////////////////////////////////////

//texi
// The @code{RedundMolecularCoor} class implements redundant sets of internal
// coordinates. @code{RedundMolecularCoor} is a
// @code{SavableState} and has @code{StateIn} and @code{KeyVal} constructors.
// @code{RedundMolecularCoor} is derived from @code{IntMolecularCoor}.
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
    //texi The KeyVal constructor (@ref{The RedundMolecularCoor KeyVal
    // Constructor}).
    RedundMolecularCoor(const RefKeyVal&);

    virtual ~RedundMolecularCoor();
    void save_data_state(StateOut&);

    //texi Actually form the variable and constant internal coordinates from
    // the simple internal coordinates.
    void form_coordinates();
    //texi Form the approximate hessian.
    void guess_hessian(RefSymmSCMatrix&hessian);
    //texi Invert the hessian.
    RefSymmSCMatrix inverse_hessian(RefSymmSCMatrix&);
};

/////////////////////////////////////////////////////////////////////////

//texi
// The @code{CartMolecularCoor} class implements Cartesian coordinates in
// a way suitable for use in geometry optimizations. @code{CartMolecularCoor}
// is a @code{SavableState} and has @code{StateIn} and @code{KeyVal}
// constructors. @code{RedundMolecularCoor} is derived from
// @code{MolecularCoor}.
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

    //texi Initializes the dimensions.
    virtual void init();
  public:
    CartMolecularCoor(RefMolecule&mol);
    CartMolecularCoor(StateIn&);
    //texi The KeyVal constructor (@ref{The CartMolecularCoor KeyVal
    // Constructor}).
    CartMolecularCoor(const RefKeyVal&);

    virtual ~CartMolecularCoor();
  
    void save_data_state(StateOut&);

    //texi These implement the virtual functions inherited from
    // @code{MolecularCoor} (@ref{The MolecularCoor Public Interface}).
    virtual RefSCDimension dim();
    virtual int to_cartesian(RefSCVector&internal);
    virtual int to_internal(RefSCVector&internal);
    virtual int to_cartesian(RefSCVector&cartesian,RefSCVector&internal);
    virtual int to_internal(RefSCVector&internal,RefSCVector&cartesian);
    virtual int to_cartesian(RefSymmSCMatrix&cart,RefSymmSCMatrix&internal);
    virtual int to_internal(RefSymmSCMatrix&internal,RefSymmSCMatrix&cart);
    virtual void print(SCostream& =SCostream::cout);
    virtual void print_simples(SCostream& =SCostream::cout);
    void guess_hessian(RefSymmSCMatrix&hessian);
    RefSymmSCMatrix inverse_hessian(RefSymmSCMatrix&);
};

#endif
