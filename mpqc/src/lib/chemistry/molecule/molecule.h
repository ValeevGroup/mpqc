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
#include <util/misc/units.h>
#include <math/symmetry/pointgrp.h>
#include <math/scmat/vector3.h>
#include <math/scmat/matrix.h>
#include <chemistry/molecule/atominfo.h>

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
    int natoms_;
    RefAtomInfo atominfo_;
    RefPointGroup pg_;
    RefUnits geometry_units_;
    double **r_;
    int *Z_;

    // these are optional
    double *mass_;
    char **labels_;

    void clear();
  public:
    Molecule();
    Molecule(Molecule&);
    Molecule(StateIn&);
    //. The \clsnmref{KeyVal} constructor.
    Molecule(const RefKeyVal&input);

    virtual ~Molecule();

    Molecule& operator=(Molecule&);

    //. Add an \clsnmref{AtomicCenter} to the \clsnm{Molecule}.
    void add_atom(int Z,double x,double y,double z,
                  const char * = 0, double mass = 0.0);

    //. Print information about the molecule.
    virtual void print(ostream& =cout);

    //. Returns the number of atoms in the molcule.
    int natom() const { return natoms_; }

    int Z(int atom) const { return Z_[atom]; }
    double &r(int atom, int xyz) { return r_[atom][xyz]; }
    const double &r(int atom, int xyz) const { return r_[atom][xyz]; }
    double *r(int atom) { return r_[atom]; }
    const double *r(int atom) const { return r_[atom]; }
    double mass(int atom) const;
    const char *label(int atom) const;

    //. Takes an (x, y, z) postion and finds an atom within the
    //given tolerance distance.  If no atom is found -1 is returned.
    int atom_at_position(double *, double tol = 0.05);

    //. Returns the index of the atom with the given \vrbl{label}.
    // If the label cannot be found \srccd{-1} is returned.
    int atom_label_to_index(const char *label) const;

    //. Returns a \srccd{double*} containing the nuclear
    //charges of the atoms.  The caller is responible for
    //freeing the return value.
    double *charges() const;

    //. Returns the total nuclear charge.
    int nuclear_charge() const;

    //. Sets the \clsnmref{PointGroup} of the molecule.
    void set_point_group(const RefPointGroup&);
    //. Returns the \clsnmref{PointGroup} of the molecule.
    const RefPointGroup point_group() const;

    //. Returns a \clsnm{SCVector3} containing the cartesian coordinates of
    // the center of mass for the molecule
    SCVector3 center_of_mass();

    //. Returns the nuclear repulsion energy for the molecule
    double nuclear_repulsion_energy();
    
    //. Compute the nuclear repulsion energy first derivative with respect
    //  to the given center. */
    void nuclear_repulsion_1der(int center, double xyz[3]);

    //. Compute the electric field due to the nuclei at the given point.
    void nuclear_efield(const double *position, double* efield);
    
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

    //. Return the number of core electrons.
    int n_core_electrons();

    //. Return the maximum atomic number.
    int max_z();

    //. Return the molecules \clsnmref{AtomInfo} object.
    RefAtomInfo atominfo() const { return atominfo_; }

    void save_data_state(StateOut&);
};
SavableState_REF_dec(Molecule);

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
