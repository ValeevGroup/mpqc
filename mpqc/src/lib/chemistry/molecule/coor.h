
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
class IntCoor: virtual public SavableState {
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
    IntCoor(const char* label = 0);
    IntCoor(const RefKeyVal&);
    IntCoor(const IntCoor&);
    IntCoor(StateIn&);
    virtual ~IntCoor();
    void save_data_state(StateOut&);
    virtual const char* label() const;
    virtual double value() const;
    virtual double preferred_value() const;
    virtual const char* ctype() const = 0; // name for coor type
#ifdef __GNUC__
    virtual void print(RefMolecule =0, SCostream& =SCostream::cout);
#else
    virtual void print();
    virtual void print(RefMolecule, SCostream& =SCostream::cout);
#endif
    virtual double force_constant(RefMolecule&) = 0;
    virtual void update_value(RefMolecule&) = 0;
    virtual void bmat(RefMolecule&,RefSCVector&bmat,double coef = 1.0) = 0;
    virtual int equivalent(RefIntCoor&) = 0;
};
SavableState_REF_dec(IntCoor);
ARRAY_dec(RefIntCoor);
SET_dec(RefIntCoor);
ARRAYSET_dec(RefIntCoor);

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
    SumIntCoor(const char *);
    SumIntCoor(const RefKeyVal&);
    SumIntCoor(StateIn&);
    ~SumIntCoor();
    void save_data_state(StateOut&);
    int n();
    void add(RefIntCoor&,double coef);
    void normalize();

    // IntCoor overrides
    double preferred_value() const;
    const char* ctype() const; // name for coor type
#ifdef __GNUC__
    void print(RefMolecule = 0, SCostream& =SCostream::cout);
#else
    void print();
    void print(RefMolecule, SCostream& =SCostream::cout);
#endif
    double force_constant(RefMolecule&);
    void update_value(RefMolecule&);
    void bmat(RefMolecule&,RefSCVector&bmat,double coef = 1.0);
    int equivalent(RefIntCoor&);
};

class SSRefSetIntCoor;
typedef class SSRefSetIntCoor RefSetIntCoor;
class SetIntCoor: virtual public SavableState {
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
    SetIntCoor(const RefKeyVal&);
    SetIntCoor(StateIn&);
    virtual ~SetIntCoor();
    void save_data_state(StateOut&);
    void add(const RefIntCoor&);
    void add(const RefSetIntCoor&);
    void clear();
    int n() const;
    RefIntCoor coor(int) const;
    virtual void fd_bmat(RefMolecule&,RefSCMatrix&);
    virtual void bmat(RefMolecule&, RefSCMatrix&);
    virtual void guess_hessian(RefMolecule&,RefSymmSCMatrix&);
#ifdef __GNUC__
    virtual void print(RefMolecule =0,SCostream& =SCostream::cout);
#else
    virtual void print();
    virtual void print(RefMolecule,SCostream& =SCostream::cout);
#endif
    virtual void update_values(RefMolecule&);
    virtual void values_to_vector(RefSCVector&);
};
SavableState_REF_dec(SetIntCoor);

class MolecularCoor: virtual public SavableState
{
#   define CLASSNAME MolecularCoor
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    RefMolecule molecule_;
  public:
    MolecularCoor(RefMolecule&);
    MolecularCoor(const RefKeyVal&);
    MolecularCoor(StateIn&);
    virtual ~MolecularCoor();
    void save_data_state(StateOut&);
    virtual void print(SCostream& =SCostream::cout) = 0;
    virtual RefSCDimension dim() = 0;
    // convert molecular coordinates to and from cartesians
    virtual void to_cartesian(RefSCVector&internal) = 0;
    virtual void to_internal(RefSCVector&internal) = 0;
    // convert the gradients
    virtual void to_cartesian(RefSCVector&cartesian,RefSCVector&internal) = 0;
    virtual void to_internal(RefSCVector&internal,RefSCVector&cartesian) = 0;
    // convert the hessian
    virtual void to_cartesian(RefSymmSCMatrix&cartesian,
                              RefSymmSCMatrix&internal) =0;
    virtual void to_internal(RefSymmSCMatrix&cartesian,
                             RefSymmSCMatrix&hessian) = 0;
    virtual void guess_hessian(RefSymmSCMatrix&hessian) = 0;
};
SavableState_REF_dec(MolecularCoor);

class IntMolecularCoor: public MolecularCoor
{
#   define CLASSNAME IntMolecularCoor
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    void form_K_matrices(RefSCDimension& dredundant,
                         RefSCDimension& dfixed,
                         RefSCMatrix& K,
                         RefSCMatrix& Kfixed,
                         int*& is_totally_symmetric);
  protected:
    RefSCDimension dim_; // corresponds to the number of variable coordinates
    RefSCDimension dnatom3_; // the number of atoms x 3
    RefSCDimension dvc_; // the number of variable + constant coordinates

    RefSetIntCoor variable_; // the variable internal coordinates
    RefSetIntCoor constant_; // the constant internal coordinates
    
    RefSetIntCoor fixed_;

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
    double scale_bonds_; // scale factor for bonds
    double scale_bends_; // scale factor for bends
    double scale_tors_;  // scale factor for tors
    double scale_outs_;  // scale factor for outs

    int nextra_bonds_;
    int* extra_bonds_;

    virtual void init();
  public:
    IntMolecularCoor(RefMolecule&mol);
    IntMolecularCoor(const RefKeyVal&);
    IntMolecularCoor(StateIn&);
    virtual ~IntMolecularCoor();
    void save_data_state(StateOut&);
    virtual void form_coordinates();
    virtual RefSCDimension dim();
    virtual void to_cartesian(RefSCVector&internal);
    virtual void to_internal(RefSCVector&internal);
    virtual void to_cartesian(RefSCVector&cartesian,RefSCVector&internal);
    virtual void to_internal(RefSCVector&internal,RefSCVector&cartesian);
    virtual void to_cartesian(RefSymmSCMatrix&cart,RefSymmSCMatrix&internal);
    virtual void to_internal(RefSymmSCMatrix&internal,RefSymmSCMatrix&cart);
    virtual void print(SCostream& =SCostream::cout);
    void guess_hessian(RefSymmSCMatrix&hessian);
};

#endif
