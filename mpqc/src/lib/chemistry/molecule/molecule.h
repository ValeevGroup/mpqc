
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

//texi
// The @code{Molecule} class provides information about the groups of atoms
// we chemists like to call molecules.  @code{Molecule} is a
// @code{SavableState} and has a @code{StateIn} constructor.  @code{Molecule}
// also has a @code{KeyVal} constructor (@ref{The Molecule KeyVal Constructor}).
class Molecule: public SavableState
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
    //texi Returns the i'th @code{AtomicCenter} in the @code{atoms} array.
    AtomicCenter& get_atom(int);
    //texi @code{const} version of the above.
    const AtomicCenter& get_atom(int) const;
  public:
    Molecule();
    Molecule(Molecule&);
    Molecule(StateIn&);
    //texi The @code{KeyVal} constructor (@ref{The Molecule KeyVal
    // Constructor}).
    Molecule(const RefKeyVal&input);

    virtual ~Molecule();

    Molecule& operator=(Molecule&);

    //texi Add an @code{AtomicCenter} to the @code{Molecule}.  The first
    // argument is the index of the atom.  You should add atoms sequentially
    // starting from zero.
    void add_atom(int,AtomicCenter&);

    //texi Print information about the molecule.
    virtual void print(SCostream& =SCostream::cout);
    virtual void print(FILE*);
    //texi Returns the number of atoms in the molcule.
    int natom() const;

    //texi Returns `1' if @code{i} is a valid @code{Pix} for this molecule.
    // Returns `0' otherwise.
    int owns(Pix i);
    //texi Returns the @code{Pix} for the first atom.
    Pix first();
    //texi Sets @code{i} to point to the next atom.  @code{i} is null if there
    // are no more atoms.
    void next(Pix& i);

    //texi Returns the @code{AtomicCenter} pointed to by @code{i}.
    AtomicCenter& operator()(Pix i);
    //texi @code{const} version of the above.
    const AtomicCenter& operator()(Pix) const;
    //texi Returns the i'th @code{AtomicCenter}.
    AtomicCenter& operator[](int i);
    //texi @code{const} version of the above.
    const AtomicCenter& operator[](int) const;
    //texi Returns the i'th @code{AtomicCenter}.
    AtomicCenter& atom(int i);
    //texi @code{const} version of the above.
    const AtomicCenter& atom(int) const;

    //texi Returns a @code{PointBag_double*} containing the nuclear charges
    // of the atoms.
    PointBag_double* charges() const;

    //texi Returns the point group of the molecule (@ref{The PointGroup Class}).
    PointGroup& point_group();
    //texi @code{const} version of the above.
    const PointGroup& point_group() const;

    //texi Returns a @code{RefPoint} containing the cartesian coordinates of
    // the center of mass for the molecule
    RefPoint center_of_mass();
    //texi If the molecule contains only symmetry unique atoms, this function
    // will generate the other, redundant atoms.
    void symmetrize();

    void save_data_state(StateOut&);
};
SavableState_REF_dec(Molecule);

/////////////////////////////////////////////////////////////////////

void mol_move_to_com(RefMolecule&);
void mol_transform_to_principal_axes(RefMolecule&, int trans_frame=1);
void mol_cleanup_molecule(RefMolecule&);

int mol_num_unique_atoms(const RefMolecule&);
int * mol_find_unique_atoms(const RefMolecule&);  // returns new'd array

#endif
