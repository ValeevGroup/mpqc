
/* simple.h -- definition of the simple internal coordinate classes
 *
 *      THIS SOFTWARE FITS THE DESCRIPTION IN THE U.S. COPYRIGHT ACT OF A
 *      "UNITED STATES GOVERNMENT WORK".  IT WAS WRITTEN AS A PART OF THE
 *      AUTHOR'S OFFICIAL DUTIES AS A GOVERNMENT EMPLOYEE.  THIS MEANS IT
 *      CANNOT BE COPYRIGHTED.  THIS SOFTWARE IS FREELY AVAILABLE TO THE
 *      PUBLIC FOR USE WITHOUT A COPYRIGHT NOTICE, AND THERE ARE NO
 *      RESTRICTIONS ON ITS USE, NOW OR SUBSEQUENTLY.
 *
 *  Author:
 *      E. T. Seidl
 *      Bldg. 12A, Rm. 2033
 *      Computer Systems Laboratory
 *      Division of Computer Research and Technology
 *      National Institutes of Health
 *      Bethesda, Maryland 20892
 *      Internet: seidl@alw.nih.gov
 *      February, 1993
 */

#ifndef _intco_simple_h
#define _intco_simple_h

#ifdef __GNUC__
#pragma interface
#endif


#include <iostream.h>
#include <stdio.h>

#include <util/class/class.h>
#include <util/state/state.h>
#include <util/keyval/keyval.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/coor.h>

#include <math/scmat/vector3.h>

//////////////////////////////////////////////////////////////////////////

//texi
// The @code{SimpleCo} class is an abstract base class for the simple
// internal coordinates (@ref{The Simple Internal Coordinate Classes}).
// @code{SimpleCo} itself inherits from @code{IntCoor} (@ref{The IntCoor
// Class}).  Like @code{IntCoor}, @code{SimpleCo} is a @code{SavableState}
// and has @code{StateIn} and @code{KeyVal} constructors.
class SimpleCo : public IntCoor {
#   define CLASSNAME SimpleCo
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    int natoms_;
    int *atoms;

  public:
    SimpleCo();
    //texi This constructor takes an integer argument which is the number of
    // atoms needed to describe the coordinate.  A second optional char*
    // argument is a label for the coordinate.  This argument is passed on
    // to the @code{IntCoor} constructor (@ref{The IntCoor Public Interface}).
    SimpleCo(int,const char* =0);
    //texi The KeyVal constructor is discussed in the next section (@ref{
    // The SimpleCo KeyVal Constructor}).
    SimpleCo(const RefKeyVal&,int);

    virtual ~SimpleCo();

    //texi Returns the number of atoms in the coordinate.
    int natoms() const;
    //texi Returns the index of the i'th atom in the coordinate.
    int operator[](int i) const;

    void save_data_state(StateOut&);
    SimpleCo(StateIn&);

    virtual int operator==(SimpleCo&);
    int operator!=(SimpleCo&u);

    // these IntCoor members are implemented in term of
    // the calc_force_con and calc_intco members.
    //texi Returns an approximate force constant (a la Almlof).
    double force_constant(RefMolecule&);
    //texi Recalculates the value of the coordinate based on the geometry
    // in the Molecule.
    void update_value(RefMolecule&);
    //texi Fill in a row of the B matrix.
    void bmat(RefMolecule&,RefSCVector&bmat,double coef = 1.0);

    //texi Calculates an approximate force constant and returns it's value.
    virtual double calc_force_con(Molecule&) = 0;
    //texi Calculate the value of the coordinate based on what's in Molecule.
    // If given a double*, fill in that part of the B matrix.  If the bmatrix
    // is to be calculated, the third argument gives the coefficient.
    virtual double calc_intco(Molecule&, double* =0, double =1) = 0;

#ifdef __GNUC__
    //texi Print the coordinate.
    void print(RefMolecule =0, SCostream& = SCostream::cout);
#else
    void print();
    void print(RefMolecule, SCostream& = SCostream::cout);
#endif
    
    //texi Tests to see if two coordinates are equivalent to each other.
    // This is false if the atoms don't match.
    int equivalent(RefIntCoor&);
  };
SavableState_REF_dec(SimpleCo);


/////////////////////////////////////////////////////////////////////////

#define SimpleCo_DECLARE(classname)					      \
  public:								      \
    virtual classname& operator=(const classname&);			      \
    SimpleCo& operator=(const SimpleCo&);				      \
    double calc_force_con(Molecule&);					      \
    double calc_intco(Molecule&, double* =0, double =1);		      \
    classname(StateIn&);						      \
    void save_data_state(StateOut&);                                          \
  private:

#define SimpleCo_IMPL_eq(classname)					      \
SimpleCo& classname::operator=(const SimpleCo& c)			      \
{									      \
  classname *cp = classname::castdown((SimpleCo*)&c);			      \
  if(cp) {								      \
      *this=*cp;							      \
    }									      \
  else {								      \
      natoms_ = 0;							      \
      atoms = 0;							      \
    }									      \
									      \
  return *this;								      \
  }

#define SimpleCo_IMPL_StateIn(classname)				      \
classname::classname(StateIn&si):					      \
  SimpleCo(si)								      \
{									      \
}

#define SimpleCo_IMPL_save_data_state(classname)			      \
void classname::save_data_state(StateOut&so)				      \
{									      \
  SimpleCo::save_data_state(so);					      \
}

#define SimpleCo_IMPL(classname)		\
        SimpleCo_IMPL_eq(classname)		\
        SimpleCo_IMPL_StateIn(classname)	\
        SimpleCo_IMPL_save_data_state(classname)

/////////////////////////////////////////////////////////////////////////

//texi
// @code{StreSimpleCo} is used to describe distances between atoms.
// 
class StreSimpleCo : public SimpleCo {
#   define CLASSNAME StreSimpleCo
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
SimpleCo_DECLARE(StreSimpleCo)
  public:
    StreSimpleCo();
    StreSimpleCo(const StreSimpleCo&);
    //texi This constructor takes a string containing a label, and two integers
    // which are the indices of the atoms we're measuring the distance between.
    // Atom numbering begins at atom 1, not atom 0.
    StreSimpleCo(const char*, int, int);
    //texi The KeyVal constructor.  This calls the @code{SimpleCo} keyval
    // constructor with an integer argument of 2 (@ref{The SimpleCo KeyVal
    // Constructor}).
    StreSimpleCo(const RefKeyVal&);

    ~StreSimpleCo();

    //texi Always returns the string "STRE".
    const char * ctype() const;

    //texi Returns the distance between the two atoms in atomic units.
    double bohr() const;
    //texi Returns the distance between the two atoms in angstrom units.
    double angstrom() const;
    //texi Returns the distance between the two atoms in angstrom units.
    double preferred_value() const;
  };

typedef StreSimpleCo Stre;

/////////////////////////////////////////////////////////////////////////

static double rtd = 180.0/3.14159265358979323846;

//texi
// @code{BendSimpleCo} is used to describe the angle abc formed by three atoms
// a, b, and c.
class BendSimpleCo : public SimpleCo { 
#   define CLASSNAME BendSimpleCo
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
SimpleCo_DECLARE(BendSimpleCo)
  public:
    BendSimpleCo();
    BendSimpleCo(const BendSimpleCo&);
    //texi This constructor takes a string containing a label, and three
    // integers a, b, and c which give the indices of the atoms involved in
    // the angle abc. Atom numbering begins at atom 1, not atom 0.
    BendSimpleCo(const char*, int, int, int);
    //texi The KeyVal constructor.  This calls the @code{SimpleCo} keyval
    // constructor with an integer argument of 3 (@ref{The SimpleCo KeyVal
    // Constructor}).
    BendSimpleCo(const RefKeyVal&);

    ~BendSimpleCo();

    //texi Always returns the string "BEND".
    const char * ctype() const;
    
    //texi Returns the value of the angle abc in radians.
    double radians() const;
    //texi Returns the value of the angle abc in degrees.
    double degrees() const;
    //texi Returns the value of the angle abc in degrees.
    double preferred_value() const;
  };

typedef BendSimpleCo Bend;

/////////////////////////////////////////////////////////////////////////

//texi
// @code{TorsSimpleCo} is used to describe the angle between the plains
// abc and bcd described by atoms a, b, c, and d.
class TorsSimpleCo : public SimpleCo { 
#   define CLASSNAME TorsSimpleCo
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
SimpleCo_DECLARE(TorsSimpleCo)
  public:
    TorsSimpleCo();
    TorsSimpleCo(const TorsSimpleCo&);
    //texi This constructor takes a string containing a label, and four
    // integers a, b, c, and d which give the indices of the atoms involved in
    // the torsion angle abcd. Atom numbering begins at atom 1, not atom 0.
    TorsSimpleCo(const char *refr, int, int, int, int);
    //texi The KeyVal constructor.  This calls the @code{SimpleCo} keyval
    // constructor with an integer argument of 4 (@ref{The SimpleCo KeyVal
    // Constructor}).
    TorsSimpleCo(const RefKeyVal&);

    ~TorsSimpleCo();

    //texi Always returns the string "TORS".
    const char * ctype() const;
    
    //texi Returns the value of the angle abc in radians.
    double radians() const;
    //texi Returns the value of the angle abc in degrees.
    double degrees() const;
    //texi Returns the value of the angle abc in degrees.
    double preferred_value() const;
  };

typedef TorsSimpleCo Tors;

/////////////////////////////////////////////////////////////////////////

//texi
// @code{ScaledTorsSimpleCo} is used to describe the angle between the plains
// abc and bcd described by atoms a, b, c, and d.  It is scaled so it makes
// sense when the abc or bcd atoms are nearly colinear.
class ScaledTorsSimpleCo : public SimpleCo { 
#   define CLASSNAME ScaledTorsSimpleCo
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
SimpleCo_DECLARE(ScaledTorsSimpleCo)
  private:
    double old_torsion_;
  public:
    ScaledTorsSimpleCo();
    ScaledTorsSimpleCo(const ScaledTorsSimpleCo&);
    //texi This constructor takes a string containing a label, and four
    // integers a, b, c, and d which give the indices of the atoms involved in
    // the torsion angle abcd. Atom numbering begins at atom 1, not atom 0.
    ScaledTorsSimpleCo(const char *refr, int, int, int, int);
    //texi The KeyVal constructor.  This calls the @code{SimpleCo} keyval
    // constructor with an integer argument of 4 (@ref{The SimpleCo KeyVal
    // Constructor}).
    ScaledTorsSimpleCo(const RefKeyVal&);

    ~ScaledTorsSimpleCo();

    //texi Always returns the string "TORS".
    const char * ctype() const;
    
    //texi Returns the value of the angle abc in radians.
    double radians() const;
    //texi Returns the value of the angle abc in degrees.
    double degrees() const;
    //texi Returns the value of the angle abc in degrees.
    double preferred_value() const;
  };

typedef ScaledTorsSimpleCo ScaledTors;

/////////////////////////////////////////////////////////////////////////

//texi
// @code{OutSimpleCo} is used to describe the out of plane angle formed by
// the bond a-b, and the plane bcd.
class OutSimpleCo : public SimpleCo { 
#   define CLASSNAME OutSimpleCo
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
SimpleCo_DECLARE(OutSimpleCo)
  public:
    OutSimpleCo();
    OutSimpleCo(const OutSimpleCo&);
    //texi This constructor takes a string containing a label, and four
    // integers a, b, c, and d which give the indices of the atoms involved in
    // the out-of-plane angle abcd. Atom numbering begins at atom 1, not
    // atom 0.
    OutSimpleCo(const char *refr, int, int, int, int);
    //texi The KeyVal constructor.  This calls the @code{SimpleCo} keyval
    // constructor with an integer argument of 4 (@ref{The SimpleCo KeyVal
    // Constructor}).
    OutSimpleCo(const RefKeyVal&);

    ~OutSimpleCo();

    //texi Always returns the string "OUT".
    const char * ctype() const;
    
    //texi Returns the value of the angle abc in radians.
    double radians() const;
    //texi Returns the value of the angle abc in degrees.
    double degrees() const;
    //texi Returns the value of the angle abc in degrees.
    double preferred_value() const;
  };

typedef OutSimpleCo Out;

/////////////////////////////////////////////////////////////////////////

//texi
// @code{LinIPSimpleCo} is used to describe the distortion of linear
// angles.  It is intended for use in finite displacement calculations.
class LinIPSimpleCo : public SimpleCo { 
#   define CLASSNAME LinIPSimpleCo
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
SimpleCo_DECLARE(LinIPSimpleCo)
  public:
    LinIPSimpleCo();
    LinIPSimpleCo(const LinIPSimpleCo&);
    //texi This constructor takes a string containing a label, and four
    // integers a, b, c, and d which give the indices of the atoms involved in
    // the linear angle abc.  d is an atom in the plane of the
    // distortion.  Atom numbering begins at atom 1, not atom 0.
    LinIPSimpleCo(const char *refr, int, int, int, int);
    //texi The KeyVal constructor.  This calls the @code{SimpleCo} keyval
    // constructor with an integer argument of 4 (@ref{The SimpleCo KeyVal
    // Constructor}).
    LinIPSimpleCo(const RefKeyVal&);

    ~LinIPSimpleCo();

    //texi Always returns the string "LINIP".
    const char * ctype() const;

    //texi Returns the value of the angle abc in radians.
    double radians() const;
    //texi Returns the value of the angle abc in degrees.
    double degrees() const;
    //texi Returns the value of the angle abc in degrees.
    double preferred_value() const;
  };

typedef LinIPSimpleCo LinIP;

/////////////////////////////////////////////////////////////////////////

//texi
// @code{LinIPSimpleCo} is used to describe the distortion of linear
// angles.  It is intended for use in finite displacement calculations.
class LinOPSimpleCo : public SimpleCo { 
#   define CLASSNAME LinOPSimpleCo
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
SimpleCo_DECLARE(LinOPSimpleCo)
  public:
    LinOPSimpleCo();
    LinOPSimpleCo(const LinOPSimpleCo&);
    //texi This constructor takes a string containing a label, and four
    // integers a, b, c, and d which give the indices of the atoms involved in
    // the linear angle abc.  d is an atom perpendicular to the plane of the
    // distortion.  Atom numbering begins at atom 1, not atom 0.
    LinOPSimpleCo(const char *refr, int =0, int =0, int =0, int =0);
    //texi The KeyVal constructor.  This calls the @code{SimpleCo} keyval
    // constructor with an integer argument of 4 (@ref{The SimpleCo KeyVal
    // Constructor}).
    LinOPSimpleCo(const RefKeyVal&);

    ~LinOPSimpleCo();

    //texi Always returns the string "LINIP".
    const char * ctype() const;

    //texi Returns the value of the angle abc in radians.
    double radians() const;
    //texi Returns the value of the angle abc in degrees.
    double degrees() const;
    //texi Returns the value of the angle abc in degrees.
    double preferred_value() const;
  };

typedef LinOPSimpleCo LinOP;

#endif /* _intco_simple_h */
