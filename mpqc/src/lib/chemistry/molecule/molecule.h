//
// molecule.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

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

    //. Takes an (x, y, z) postion and finds an atom within the
    //given tolerance distance.  If no atom is found -1 is returned.
    int atom_at_position(double *, double tol = 0.05);

    //. Returns the index of the atom with the given \vrbl{label}.
    // If the label cannot be found \srccd{-1} is returned.
    int atom_label_to_index(const char *label) const;

    //. Returns a \srccd{\clsnm{PointBag\_double}*} containing the nuclear
    //charges of the atoms.
    PointBag_double* charges() const;

    //. Sets the \clsnmref{PointGroup} of the molecule.
    void set_point_group(const PointGroup&);
    //. Returns the \clsnmref{PointGroup} of the molecule.
    const PointGroup& point_group() const;

    //. Returns a \clsnm{RefPoint} containing the cartesian coordinates of
    // the center of mass for the molecule
    RefPoint center_of_mass();

    //. Returns the nuclear repulsion energy for the molecule
    double nuclear_repulsion_energy();
    
    //. Compute the nuclear repulsion energy first derivative with respect
    //  to the given center. */
    void nuclear_repulsion_1der(int center, double xyz[3]);
    
    //. If the molecule contains only symmetry unique atoms, this function
    // will generate the other, redundant atoms.
    void symmetrize();

    void move_to_com();
    void transform_to_principal_axes(int trans_frame=1);
    void cleanup_molecule();
    void print_pdb(ostream& =cout, char *title =0);

    //. Compute the principal moments of inertia and, possibly, the
    // principal axes.
    void principal_moments_of_inertia(double *evals, double **evecs=0);

    int num_unique_atoms();
    int *find_unique_atoms();  // returns new'd array

    void save_data_state(StateOut&);
};
SavableState_REF_dec(Molecule);

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
