
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

#undef V_BASE
#define V_BASE virtual public SavableState

//////////////////////////////////////////////////////////////////////////

class SimpleCo : public IntCoor {
#   define CLASSNAME SimpleCo
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    int natoms_;
    int *atoms;

  public:
    SimpleCo();
    SimpleCo(int,const char* =0);
    SimpleCo(const RefKeyVal&,int);
    virtual ~SimpleCo();

    int natoms() const;
    int operator[](int i) const;

    void save_data_state(StateOut&);
    SimpleCo(StateIn&);

    virtual int operator==(SimpleCo&);
    int operator!=(SimpleCo&u);

    // these IntCoor members are implemented in term of
    // the calc_force_con and calc_intco members.
    double force_constant(RefMolecule&);
    void update_value(RefMolecule&);
    void bmat(RefMolecule&,RefSCVector&bmat,double coef = 1.0);

    virtual double calc_force_con(Molecule&) = 0;
    virtual double calc_intco(Molecule&, double* =0, double =1) = 0;

#ifdef __GNUC__
    void print(RefMolecule =0, SCostream& = SCostream::cout);
#else
    void print();
    void print(RefMolecule, SCostream& = SCostream::cout);
#endif
    
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

#define SimpleCo_IMPL(classname)					      \
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
  }									      \
classname::classname(StateIn&si):					      \
  SimpleCo(si)								      \
{									      \
}									      \
void classname::save_data_state(StateOut&so)				      \
{									      \
  SimpleCo::save_data_state(so);					      \
}

/////////////////////////////////////////////////////////////////////////

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
    StreSimpleCo(const char*, int, int);
    StreSimpleCo(const RefKeyVal&);
    //StreSimpleCo(KeyVal*,const char*,int=0);
    ~StreSimpleCo();

    const char * ctype() const;
    
    double bohr() const;
    double angstrom() const;
    double preferred_value() const;
  };

typedef StreSimpleCo Stre;

/////////////////////////////////////////////////////////////////////////

static double rtd = 180.0/3.14159265358979323846;

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
    BendSimpleCo(const char*, int, int, int);
    BendSimpleCo(const RefKeyVal&);
    //BendSimpleCo(KeyVal*,const char*,int=0);
    ~BendSimpleCo();

    const char * ctype() const;
    
    double radians() const;
    double degrees() const;
    double preferred_value() const;
  };

typedef BendSimpleCo Bend;

/////////////////////////////////////////////////////////////////////////

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
    TorsSimpleCo(const char *refr, int, int, int, int);
    TorsSimpleCo(const RefKeyVal&);
    //TorsSimpleCo(KeyVal*,const char*,int=0);
    ~TorsSimpleCo();

    const char * ctype() const;
    
    double radians() const;
    double degrees() const;
    double preferred_value() const;
  };

typedef TorsSimpleCo Tors;

/////////////////////////////////////////////////////////////////////////

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
    OutSimpleCo(const char *refr, int, int, int, int);
    OutSimpleCo(const RefKeyVal&);
    //OutSimpleCo(KeyVal*,const char*,int=0);
    ~OutSimpleCo();

    const char * ctype() const;
    
    double radians() const;
    double degrees() const;
    double preferred_value() const;
  };

typedef OutSimpleCo Out;

/////////////////////////////////////////////////////////////////////////

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
    LinIPSimpleCo(const char *refr, int, int, int, int);
    LinIPSimpleCo(const RefKeyVal&);
    //LinIPSimpleCo(KeyVal*,const char*,int=0);
    ~LinIPSimpleCo();

    const char * ctype() const;

    void set_theta(const double*, const double*, const double*, const double*);

    double radians() const;
    double degrees() const;
    double preferred_value() const;
  };

typedef LinIPSimpleCo LinIP;

/////////////////////////////////////////////////////////////////////////

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
    LinOPSimpleCo(const char *refr, int =0, int =0, int =0, int =0);
    LinOPSimpleCo(const RefKeyVal&);
    //LinOPSimpleCo(KeyVal*,const char*,int=0);
    ~LinOPSimpleCo();

    const char * ctype() const;

    void set_theta(const double*, const double*, const double*, const double*);

    double radians() const;
    double degrees() const;
    double preferred_value() const;
  };

typedef LinOPSimpleCo LinOP;

/////////////////////////////////////////////////////////////////////

/*
 * these are some utility routines
 */

/////////////////////////////////////////////////////////////////////

#undef V_BASE

#endif /* _intco_simple_h */
