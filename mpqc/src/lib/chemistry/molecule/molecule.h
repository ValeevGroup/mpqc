
#ifndef _chemistry_molecule_molecule_h
#define _chemistry_molecule_molecule_h

#ifdef __GNUC__
#pragma interface
#endif

#include <stdio.h>
#include <Pix.h>
#ifdef __GNUC__
#include <ostream.h>
#endif
#include <util/misc/scostream.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/keyval/keyval.h>
#include <math/topology/point.h>
#include <math/topology/pointbag.h>
#include <math/symmetry/pointgrp.h>
#include <math/scmat/matrix.h>
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
    PointGroup pg;
    AtomicCenter* atoms;
    int natoms;
    AtomicCenter& get_atom(int);
    const AtomicCenter& get_atom(int) const;
    RefSCDimension dnatom3_; // number of atoms times 3
  public:
    Molecule();
    Molecule(Molecule&);
    Molecule(StateIn&);
    Molecule(KeyVal&input);
    virtual ~Molecule();

    RefSCDimension dim_natom3(); // return natom3_;

    Molecule& operator=(Molecule&);

    void add_atom(int,AtomicCenter&);

    virtual void print(SCostream& =SCostream::cout);
    virtual void print(FILE*);
    int natom() const;

    int owns(Pix);
    Pix first();
    void next(Pix&);

    AtomicCenter& operator()(Pix);
    AtomicCenter& operator[](int);
    AtomicCenter& atom(int);
    const AtomicCenter& operator()(Pix) const;
    const AtomicCenter& operator[](int) const;
    const AtomicCenter& atom(int) const;

    PointBag_double* charges() const;

    PointGroup& point_group();
    const PointGroup& point_group() const;
    Point center_of_mass();

    void save_data_state(StateOut&);
};
SavableState_REF_dec(Molecule);

/////////////////////////////////////////////////////////////////////

// this may or may not increase the number of atoms
void mol_symmetrize_molecule(Molecule&);

// these change the coordinates of the molecule
void mol_move_to_com(Molecule&);
void mol_transform_to_principal_axes(Molecule&);
void mol_cleanup_molecule(Molecule&);

// these have no side effects
int mol_num_unique_atoms(Molecule&);
int * mol_find_unique_atoms(Molecule&);  // returns new'd array

#endif
