
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

//.  The \clsnm{Molecule} class provides information about the groups of
//atoms we chemists like to call molecules.  \clsnm{Molecule} is a
//\clsnmref{SavableState} and has a \clsnmref{StateIn} constructor.
//\clsnm{Molecule} also has a \clsnmref{KeyVal} constructor.
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
    //. Returns the \vrbl{i}'th \clsnmref{AtomicCenter}.
    AtomicCenter& get_atom(int i);
    //. \srccd{const} version of the above.
    const AtomicCenter& get_atom(int) const;
  public:
    Molecule();
    Molecule(Molecule&);
    Molecule(StateIn&);
    //. The \clsnmref{KeyVal} constructor.
    Molecule(const RefKeyVal&input);

    virtual ~Molecule();

    Molecule& operator=(Molecule&);

    //. Add an \clsnmref{AtomicCenter} to the \clsnm{Molecule}.  The first
    // argument is the index of the atom.  You should add atoms sequentially
    // starting from zero.
    void add_atom(int,AtomicCenter&);

    //. Print information about the molecule.
    virtual void print(ostream& =cout);
    virtual void print(FILE*);
    //. Returns the number of atoms in the molcule.
    int natom() const;

    //. Returns `1' if \vrbl{i} is a valid \clsnm{Pix} for this molecule.
    // Returns `0' otherwise.
    int owns(Pix i);
    //. Returns the \clsnm{Pix} for the first atom.
    Pix first();
    //. Sets \vrbl{i} to point to the next atom.  \vrbl{i} is null if there
    // are no more atoms.
    void next(Pix& i);

    //. Returns the \clsnmref{AtomicCenter} pointed to by \vrbl{i}.
    AtomicCenter& operator()(Pix i);
    //. \srccd{const} version of the above.
    const AtomicCenter& operator()(Pix) const;
    //. Returns the i'th \clsnmref{AtomicCenter}.
    AtomicCenter& operator[](int i);
    //. \srccd{const} version of the above.
    const AtomicCenter& operator[](int) const;
    //. Returns the i'th \clsnmref{AtomicCenter}.
    AtomicCenter& atom(int i);
    //. \srccd{const} version of the above.
    const AtomicCenter& atom(int) const;

    //. Returns the index of the atom with the given \vrbl{label}.
    // If the label cannot be found \srccd{-1} is returned.
    int atom_label_to_index(const char *label) const;

    //. Returns a \srccd{\clsnm{PointBag\_double}*} containing the nuclear
    //charges of the atoms.
    PointBag_double* charges() const;

    //. Returns the \clsnmref{PointGroup} of the molecule.
    PointGroup& point_group();
    //. \srccd{const} version of the above.
    const PointGroup& point_group() const;

    //. Returns a \clsnm{RefPoint} containing the cartesian coordinates of
    // the center of mass for the molecule
    RefPoint center_of_mass();

    //. Returns the nuclear repulsion energy for the molecule
    double nuclear_repulsion_energy();
    
    //. If the molecule contains only symmetry unique atoms, this function
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
