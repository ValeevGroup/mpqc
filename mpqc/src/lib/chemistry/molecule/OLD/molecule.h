
#ifndef _molecule_h
#define _molecule_h

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
  virtual public DescribedClass, virtual public SavableState
{
  DescribedClass_DECLARE(Molecule)
  SavableState_DECLARE(Molecule)
  private:
  //AtomicCenterXPlex atoms;
    AtomicCenter* atoms;
    int natoms;
    AtomicCenter& get_atom(int);
    const AtomicCenter& get_atom(int) const;
  public:
    Molecule();
    Molecule(Molecule&);
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
    void restore_data_state(int,StateIn&);
};

#endif
