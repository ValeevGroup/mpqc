
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


#include <iostream.h>
#include <stdio.h>

#include <util/class/class.h>
#include <util/state/state.h>
#include <util/keyval/keyval.h>
#include <chemistry/molecule/molecule.h>

#undef V_BASE
#define V_BASE virtual public DescribedClass, virtual public SavableState

//////////////////////////////////////////////////////////////////////////

class SimpleCoPtr;
class SimpleCoList;

class SCCount {
  friend class SimpleCoPtr;

  private:
    int count;
  protected:
    SCCount(): count(0) {}
  };

//////////////////////////////////////////////////////////////////////////

class SimpleCo : public SCCount, V_BASE {
DescribedClass_DECLARE_ABSTRACT(SimpleCo)
SavableState_DECLARE(SimpleCo)
  protected:
    char *ref;
    int natoms_;
    int *atoms;
    double value_;

    SimpleCo();
    SimpleCo(int,const char* =0);

  public:
    virtual ~SimpleCo();

    inline const char * reference() const { return ref; }
    inline const int natoms() const { return natoms_; }
    inline const int operator[](int i) const { return atoms[i]; }
    inline double value() const { return value_; }
    virtual double preferred_value() const =0;

    virtual const char * ctype() const =0;

    virtual void init();

    void save_data_state(StateOut&);
    void restore_data_state(int,StateIn&);

    virtual int operator==(SimpleCo&);
    inline int operator!=(SimpleCo&u) { return !(*this==u); }

    virtual double calc_force_con(Molecule&) =0;
    virtual double calc_intco(Molecule&, double* =0, double =1) =0;

    virtual void print(ostream&, const char* =" ") const =0;
    virtual void print(FILE* =stdout, const char* =" ") const =0;

    friend ostream& operator<<(ostream&,SimpleCo&);
  };

/////////////////////////////////////////////////////////////////////////

#define SimpleCo_DECLARE(classname) \
  public: \
    virtual classname& operator=(const classname&); \
    SimpleCo& operator=(const SimpleCo&); \
    double calc_force_con(Molecule&); \
    double calc_intco(Molecule&, double* =0, double =1); \
    void print(ostream&, const char* =" ") const; \
    void print(FILE* =stdout, const char* =" ") const; \
    inline void save_data_state(StateOut&so) \
                             { SimpleCo::save_parent_state(so); } \
    inline void restore_data_state(int,StateIn&si) \
                             { SimpleCo::restore_parent_state(si); } \
  private:

#define SimpleCo_IMPL(classname) \
SimpleCo& classname::operator=(const SimpleCo& c) \
{ \
  classname *cp = classname::castdown((SimpleCo*)&c); \
  if(cp) \
    *this=*cp; \
  else \
    init(); \
 \
  return *this; \
  }

/////////////////////////////////////////////////////////////////////////

class StreSimpleCo : public SimpleCo, V_BASE {
DescribedClass_DECLARE(StreSimpleCo)
SavableState_DECLARE(StreSimpleCo)
SimpleCo_DECLARE(StreSimpleCo)
  public:
    StreSimpleCo();
    StreSimpleCo(const StreSimpleCo&);
    StreSimpleCo(const char*, int, int);
    StreSimpleCo(KeyVal*,const char*,int=0);
    ~StreSimpleCo();

    inline const char * ctype() const { return "STRE"; }

    inline double bohr() const { return value_; }
    inline double angstrom() const { return value_*0.52917706; }
    inline double preferred_value() const { return value_*0.52917706; }
  };

typedef StreSimpleCo Stre;

/////////////////////////////////////////////////////////////////////////

static double rtd = 180.0/3.14159265358979323846;

class BendSimpleCo : public SimpleCo, V_BASE { 
DescribedClass_DECLARE(BendSimpleCo)
SavableState_DECLARE(BendSimpleCo)
SimpleCo_DECLARE(BendSimpleCo)
  public:
    BendSimpleCo();
    BendSimpleCo(const BendSimpleCo&);
    BendSimpleCo(const char*, int, int, int);
    BendSimpleCo(KeyVal*,const char*,int=0);
    ~BendSimpleCo();

    inline const char * ctype() const { return "BEND"; }

    inline double radians() const { return value_; }
    inline double degrees() const { return value_*rtd; }
    inline double preferred_value() const { return value_*rtd; }
  };

typedef BendSimpleCo Bend;

/////////////////////////////////////////////////////////////////////////

class TorsSimpleCo : public SimpleCo, V_BASE { 
DescribedClass_DECLARE(TorsSimpleCo)
SavableState_DECLARE(TorsSimpleCo)
SimpleCo_DECLARE(TorsSimpleCo)
  public:
    TorsSimpleCo();
    TorsSimpleCo(const TorsSimpleCo&);
    TorsSimpleCo(const char *refr, int, int, int, int);
    TorsSimpleCo(KeyVal*,const char*,int=0);
    ~TorsSimpleCo();

    inline const char * ctype() const { return "TORS"; }

    inline double radians() const { return value_; }
    inline double degrees() const { return value_*rtd; }
    inline double preferred_value() const { return value_*rtd; }
  };

typedef TorsSimpleCo Tors;

/////////////////////////////////////////////////////////////////////////

class OutSimpleCo : public SimpleCo, V_BASE { 
DescribedClass_DECLARE(OutSimpleCo)
SavableState_DECLARE(OutSimpleCo)
SimpleCo_DECLARE(OutSimpleCo)
  public:
    OutSimpleCo();
    OutSimpleCo(const OutSimpleCo&);
    OutSimpleCo(const char *refr, int, int, int, int);
    OutSimpleCo(KeyVal*,const char*,int=0);
    ~OutSimpleCo();

    inline const char * ctype() const { return "OUT"; }

    inline double radians() const { return value_; }
    inline double degrees() const { return value_*rtd; }
    inline double preferred_value() const { return value_*rtd; }
  };

typedef OutSimpleCo Out;

/////////////////////////////////////////////////////////////////////////

class LinIPSimpleCo : public SimpleCo, V_BASE { 
DescribedClass_DECLARE(LinIPSimpleCo)
SavableState_DECLARE(LinIPSimpleCo)
SimpleCo_DECLARE(LinIPSimpleCo)
  public:
    LinIPSimpleCo();
    LinIPSimpleCo(const LinIPSimpleCo&);
    LinIPSimpleCo(const char *refr, int, int, int, int);
    LinIPSimpleCo(KeyVal*,const char*,int=0);
    ~LinIPSimpleCo();

    inline const char * ctype() const { return "LINIP"; }

    void set_theta(const double*, const double*, const double*, const double*);

    inline double radians() const { return value_; }
    inline double degrees() const { return value_*rtd; }
    inline double preferred_value() const { return value_*rtd; }
  };

typedef LinIPSimpleCo LinIP;

/////////////////////////////////////////////////////////////////////////

class LinOPSimpleCo : public SimpleCo, V_BASE { 
DescribedClass_DECLARE(LinOPSimpleCo)
SavableState_DECLARE(LinOPSimpleCo)
SimpleCo_DECLARE(LinOPSimpleCo)
  public:
    LinOPSimpleCo();
    LinOPSimpleCo(const LinOPSimpleCo&);
    LinOPSimpleCo(const char *refr, int =0, int =0, int =0, int =0);
    LinOPSimpleCo(KeyVal*,const char*,int=0);
    ~LinOPSimpleCo();

    inline const char * ctype() const { return "LINOP"; }

    void set_theta(const double*, const double*, const double*, const double*);

    inline double radians() const { return value_; }
    inline double degrees() const { return value_*rtd; }
    inline double preferred_value() const { return value_*rtd; }
  };

typedef LinOPSimpleCo LinOP;

/////////////////////////////////////////////////////////////////////

/*
 * these are some utility routines
 */

SimpleCoList * Geom_read_simples(KeyVal*);
SimpleCoList * Geom_form_simples(Molecule&);
void Geom_calc_simples(SimpleCoList*,Molecule&);

void Geom_print_pretty(SimpleCoList*);
void Geom_print_pretty(ostream&,SimpleCoList*,const double* =0);

/////////////////////////////////////////////////////////////////////

#undef V_BASE

#endif /* _intco_simple_h */
