
#ifndef _chemistry_molecule_molecule_h
#define _chemistry_molecule_molecule_h

#include <stdio.h>
#include <Pix.h>

#include <util/class/class.h>
#include <util/state/state.h>
#include <util/keyval/keyval.h>
#include <math/topology/point.h>
#include <math/topology/pointbag.h>
#include <chemistry/molecule/chemelem.h>
#include <chemistry/molecule/atomcent.h>

class PointBag_double;

//#include "atomcentXPlex.h"

// Generic Molecule
class Molecule :
  virtual public SavableState
{
#   define CLASSNAME Molecule
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
  //AtomicCenterXPlex atoms;
    AtomicCenter* atoms;
    int natoms;
    AtomicCenter& get_atom(int);
    const AtomicCenter& get_atom(int) const;
  public:
    Molecule();
    Molecule(Molecule&);
    Molecule(StateIn&);
    Molecule(KeyVal&input);
    virtual ~Molecule();

    Molecule& operator=(Molecule&);

    void add_atom(int,AtomicCenter&);

    virtual void print(FILE*fp=stdout);
    int natom() const;
    int owns(Pix);
    Pix first();
    void next(Pix&);
    AtomicCenter& operator()(Pix);
    AtomicCenter& operator[](int);
    const AtomicCenter& operator()(Pix) const; 
    const AtomicCenter& operator[](int) const;
    PointBag_double* charges() const;

    void save_data_state(StateOut&);
};
DescribedClass_REF_dec(Molecule);

#endif
