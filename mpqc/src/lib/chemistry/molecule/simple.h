
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

#include <util/class/class.h>
#include <util/state/state.h>
#include <util/keyval/keyval.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/coor.h>

#include <math/scmat/vector3.h>

//////////////////////////////////////////////////////////////////////////

//.  The \clsnm{SimpleCo} class is an abstract base class for the simple
//internal coordinates.  \clsnm{SimpleCo} itself inherits from
//\clsnmref{IntCoor}.  Like \clsnmref{IntCoor}, \clsnm{SimpleCo} is a
//\clsnmref{SavableState} and has \clsnmref{StateIn} and \clsnmref{KeyVal}
//constructors.
class SimpleCo : public IntCoor {
#   define CLASSNAME SimpleCo
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    int natoms_;
    int *atoms;

  public:
    SimpleCo();
    //. This constructor takes an integer argument which is the number of
    // atoms needed to describe the coordinate.  A second optional char*
    // argument is a label for the coordinate.  This argument is passed on
    // to the \clsnmref{IntCoor} constructor.
    SimpleCo(int,const char* =0);
    //. The \clsnmref{KeyVal} constructor requires the number of
    // atoms.
    SimpleCo(const RefKeyVal&,int natom);

    virtual ~SimpleCo();

    //. Returns the number of atoms in the coordinate.
    int natoms() const;
    //. Returns the index of the i'th atom in the coordinate.
    int operator[](int i) const;

    void save_data_state(StateOut&);
    SimpleCo(StateIn&);

    virtual int operator==(SimpleCo&);
    int operator!=(SimpleCo&u);

    // these IntCoor members are implemented in term of
    // the calc_force_con and calc_intco members.
    //. Returns an approximate force constant (a la Almlof).
    double force_constant(RefMolecule&);
    //. Recalculates the value of the coordinate based on the geometry
    // in the Molecule.
    void update_value(const RefMolecule&);
    //. Fill in a row of the B matrix.
    void bmat(const RefMolecule&,RefSCVector&bmat,double coef = 1.0);

    //. Calculates an approximate force constant and returns it's value.
    virtual double calc_force_con(Molecule&) = 0;
    //. Calculate the value of the coordinate based on what's in Molecule.
    // If given a double*, fill in that part of the B matrix.  If the bmatrix
    // is to be calculated, the third argument gives the coefficient.
    virtual double calc_intco(Molecule&, double* =0, double =1) = 0;

#ifdef __GNUC__
    //. Print the coordinate.
    void print(RefMolecule =0, ostream& = cout);
#else
    void print();
    void print(RefMolecule, ostream& = cout);
#endif
    
    //. Tests to see if two coordinates are equivalent to each other.
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

//.
// \clsnm{StreSimpleCo} is used to describe distances between atoms.
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
    //. This constructor takes a string containing a label, and two integers
    // which are the indices of the atoms we're measuring the distance between.
    // Atom numbering begins at atom 1, not atom 0.
    StreSimpleCo(const char*, int, int);
    //. The \clsnmref{KeyVal} constructor.  This calls the
    //\clsnmref{SimpleCo} keyval constructor with an integer argument of 2.
    StreSimpleCo(const RefKeyVal&);

    ~StreSimpleCo();

    //. Always returns the string "STRE".
    const char * ctype() const;

    //. Returns the distance between the two atoms in atomic units.
    double bohr() const;
    //. Returns the distance between the two atoms in angstrom units.
    double angstrom() const;
    //. Returns the distance between the two atoms in angstrom units.
    double preferred_value() const;
  };

typedef StreSimpleCo Stre;

/////////////////////////////////////////////////////////////////////////

static const double rtd = 180.0/3.14159265358979323846;

//.
// \clsnm{BendSimpleCo} is used to describe the angle abc formed by three atoms
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
    //. This constructor takes a string containing a label, and three
    // integers a, b, and c which give the indices of the atoms involved in
    // the angle abc. Atom numbering begins at atom 1, not atom 0.
    BendSimpleCo(const char*, int, int, int);
    //. The \clsnmref{KeyVal} constructor.  This calls the
    //\clsnmref{SimpleCo} keyval constructor with an integer argument of 3.
    BendSimpleCo(const RefKeyVal&);

    ~BendSimpleCo();

    //. Always returns the string "BEND".
    const char * ctype() const;
    
    //. Returns the value of the angle abc in radians.
    double radians() const;
    //. Returns the value of the angle abc in degrees.
    double degrees() const;
    //. Returns the value of the angle abc in degrees.
    double preferred_value() const;
  };

typedef BendSimpleCo Bend;

/////////////////////////////////////////////////////////////////////////

//.
// \clsnm{TorsSimpleCo} is used to describe the angle between the plains
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
    //. This constructor takes a string containing a label, and four
    // integers a, b, c, and d which give the indices of the atoms involved in
    // the torsion angle abcd. Atom numbering begins at atom 1, not atom 0.
    TorsSimpleCo(const char *refr, int, int, int, int);
    //. The \clsnmref{KeyVal} constructor.  This calls the
    //\clsnmref{SimpleCo} keyval constructor with an integer argument of 4.
    TorsSimpleCo(const RefKeyVal&);

    ~TorsSimpleCo();

    //. Always returns the string "TORS".
    const char * ctype() const;
    
    //. Returns the value of the angle abc in radians.
    double radians() const;
    //. Returns the value of the angle abc in degrees.
    double degrees() const;
    //. Returns the value of the angle abc in degrees.
    double preferred_value() const;
  };

typedef TorsSimpleCo Tors;

/////////////////////////////////////////////////////////////////////////

//.
// \clsnm{ScaledTorsSimpleCo} is used to describe the angle between the plains
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
    //. This constructor takes a string containing a label, and four
    // integers a, b, c, and d which give the indices of the atoms involved in
    // the torsion angle abcd. Atom numbering begins at atom 1, not atom 0.
    ScaledTorsSimpleCo(const char *refr, int, int, int, int);
    //. The \clsnmref{KeyVal} constructor.  This calls the
    //\clsnmref{SimpleCo} keyval constructor with an integer argument of 4.
    ScaledTorsSimpleCo(const RefKeyVal&);

    ~ScaledTorsSimpleCo();

    //. Always returns the string "TORS".
    const char * ctype() const;
    
    //. Returns the value of the angle abc in radians.
    double radians() const;
    //. Returns the value of the angle abc in degrees.
    double degrees() const;
    //. Returns the value of the angle abc in degrees.
    double preferred_value() const;
  };

typedef ScaledTorsSimpleCo ScaledTors;

/////////////////////////////////////////////////////////////////////////

//.
// \clsnm{OutSimpleCo} is used to describe the out of plane angle formed by
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
    //. This constructor takes a string containing a label, and four
    // integers a, b, c, and d which give the indices of the atoms involved in
    // the out-of-plane angle abcd. Atom numbering begins at atom 1, not
    // atom 0.
    OutSimpleCo(const char *refr, int, int, int, int);
    //. The \clsnmref{KeyVal} constructor.  This calls the \clsnmref{SimpleCo} keyval
    // constructor with an integer argument of 4.
    OutSimpleCo(const RefKeyVal&);

    ~OutSimpleCo();

    //. Always returns the string "OUT".
    const char * ctype() const;
    
    //. Returns the value of the angle abc in radians.
    double radians() const;
    //. Returns the value of the angle abc in degrees.
    double degrees() const;
    //. Returns the value of the angle abc in degrees.
    double preferred_value() const;
  };

typedef OutSimpleCo Out;

/////////////////////////////////////////////////////////////////////////

//.
// \clsnm{LinIPSimpleCo} is used to describe the distortion of linear
// angles.
class LinIPSimpleCo : public SimpleCo { 
#   define CLASSNAME LinIPSimpleCo
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
SimpleCo_DECLARE(LinIPSimpleCo)
  private:
    SCVector3 u2;
  public:
    LinIPSimpleCo();
    LinIPSimpleCo(const LinIPSimpleCo&);
    //. This constructor takes a string containing a label, and three
    // integers a, b, and d which give the indices of the atoms involved in
    // the linear angle abc.  The last argument, u, is a unit vector
    // used to defined the direction in which distortion is measured.
    // Atom numbering begins at atom 1, not atom 0.
    LinIPSimpleCo(const char *refr, int, int, int, const SCVector3 &u);
    //. The \clsnmref{KeyVal} constructor.  This calls the \clsnm{SimpleCo} keyval
    // constructor with an integer argument of 3.
    LinIPSimpleCo(const RefKeyVal&);

    ~LinIPSimpleCo();

    //. Always returns the string "LINIP".
    const char * ctype() const;

    //. Returns the value of the angle abc in radians.
    double radians() const;
    //. Returns the value of the angle abc in degrees.
    double degrees() const;
    //. Returns the value of the angle abc in degrees.
    double preferred_value() const;
  };

typedef LinIPSimpleCo LinIP;

/////////////////////////////////////////////////////////////////////////

//.
// \clsnm{LinOPSimpleCo} is used to describe the distortion of linear
// angles.
class LinOPSimpleCo : public SimpleCo { 
#   define CLASSNAME LinOPSimpleCo
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
SimpleCo_DECLARE(LinOPSimpleCo)
  private:
    SCVector3 u2;
  public:
    LinOPSimpleCo();
    LinOPSimpleCo(const LinOPSimpleCo&);
    //. This constructor takes a string containing a label, and three
    // integers a, b, and c which give the indices of the atoms involved in
    // the linear angle abc.  The last argument, u, is a unit vector used to
    // defined the direction perpendicular to the direction in which
    // distortion is measured.  Atom numbering begins at atom 1, not atom 0.
    LinOPSimpleCo(const char *refr, int, int, int, const SCVector3 &u);
    //. The \clsnmref{KeyVal} constructor.  This calls the
    //\clsnmref{SimpleCo} keyval constructor with an integer argument of 3.
    LinOPSimpleCo(const RefKeyVal&);

    ~LinOPSimpleCo();

    //. Always returns the string "LINIP".
    const char * ctype() const;

    //. Returns the value of the angle abc in radians.
    double radians() const;
    //. Returns the value of the angle abc in degrees.
    double degrees() const;
    //. Returns the value of the angle abc in degrees.
    double preferred_value() const;
  };

typedef LinOPSimpleCo LinOP;

#endif /* _intco_simple_h */

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
