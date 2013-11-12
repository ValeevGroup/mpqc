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

#include <stdio.h>
#include <iostream>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/keyval/keyval.h>
#include <util/misc/units.h>
#include <math/symmetry/pointgrp.h>
#include <math/scmat/vector3.h>
#include <math/scmat/matrix.h>
#include <chemistry/molecule/atominfo.h>
#include <chemistry/molecule/atom.h>
#include <util/misc/xml.h>

using boost::property_tree::ptree;
namespace sc {

  /// @addtogroup ChemistryMolecule
  /// @{

/**
The Molecule class contains information about molecules.  It has a
KeyVal constructor that can create a new molecule from either a
PDB file or from a list of Cartesian coordinates.

The following ParsedKeyVal input reads from the PDB
file <tt>h2o.pdb</tt>:
<pre>
molecule<Molecule>: (
   pdb_file = "h2o.pdb"
 )
</pre>

The following input explicitly gives the atom coordinates, using the
ParsedKeyVal table notation:
<pre>
molecule<Molecule>: (
    unit=angstrom
    { atom_labels atoms           geometry            } = {
          O1         O   [ 0.000000000 0  0.369372944 ]
          H1         H   [ 0.783975899 0 -0.184686472 ]
          H2         H   [-0.783975899 0 -0.184686472 ]
     }
    )
  )
</pre>
The default units are Bohr which can be overridden with
<tt>unit=angstrom</tt>.  The <tt>atom_labels</tt> array can be
omitted.  The <tt>atoms</tt> and <tt>geometry</tt> arrays
are required.

As a special case, an atom can be given with the symbol <tt>Q</tt> or the
name <tt>charge</tt>.  Such centers are treated as point charges and not
given basis functions.  The values of the charges must be specified with a
<tt>charge</tt> vector in the Molecule input.  Since the charge vector
assign charges to all centers, including atoms, it is easiest to place all
point charge centers first in the geometry, and then give a charge vector
with a number of elements equal to the number of point charges.  The
following example shows a water molecule interacting with a point charge
having value 0.1:
<pre>
molecule<Molecule>: (
    unit=angstrom
    charge = [ 0.1 ]
    { atom_labels atoms           geometry            } = {
          Q1         Q   [ 0.0         0 10.0         ]
          O1         O   [ 0.000000000 0  0.369372944 ]
          H1         H   [ 0.783975899 0 -0.184686472 ]
          H2         H   [-0.783975899 0 -0.184686472 ]
     }
    )
  )
</pre>

This feature is designed for doing QM/MM calculations, so, by default,
methods will not include interactions between the <tt>Q</tt> centers when
computing the energy or the gradient.  To include these interactions, set
<tt>include_qq=1</tt>.

The Molecule class has a PointGroup
member object, which also has a KeyVal constructor
that is called when a Molecule is made.  The
following example constructs a molecule with \f$C_{2v}\f$ symmetry:
<pre>
molecule<Molecule>: (
    symmetry=c2v
    unit=angstrom
    { atoms         geometry            } = {
        O   [0.000000000 0  0.369372944 ]
        H   [0.783975899 0 -0.184686472 ]
     }
    )
  )
</pre>
Only the symmetry unique atoms need to be specified.  Nonunique
atoms can be given too, however, numerical errors in the
geometry specification can result in the generation of extra
atoms so be careful.
*/
class Molecule: public SavableState, virtual public DescribedXMLWritable
{
  protected:
    std::vector<Atom> atoms_;
    Ref<AtomInfo> atominfo_;
    Ref<PointGroup> pg_;
    Ref<Units> geometry_units_;
    double ref_origin_[3];   //< position of the origin of the reference coordinate system

    // symmetry equiv info
    int nuniq_;
    int *nequiv_;
    int **equiv_;
    int *atom_to_uniq_;
    void init_symmetry_info(double tol=0.5);
    void clear_symmetry_info();

    // The Z that represents a "Q" type atom.
    int q_Z_;

    // If true, include the q terms in the charge and efield routines
    bool include_q_;

    // If true, include the coupling between q-q pairs when
    // computing nuclear repulsion energy and gradients.
    bool include_qq_;

    // These vectors contain the atom indices of atoms that are not type
    // "Q" and those that are.
    std::vector<int> q_atoms_;
    std::vector<int> non_q_atoms_;

    void clear();

    // Throw an exception if an atom is duplicated.  The
    // atoms in the range [begin, natom_) are checked.
    void throw_if_atom_duplicated(int begin=0, double tol = 1e-3);
  public:
    Molecule();
    Molecule(const Molecule&);
    Molecule(StateIn&);
    /** The Molecule KeyVal constructor is used to generate a Molecule
        object from the input.  Several examples are given in the Molecule
        class overview.  The full list of keywords that are accepted is
        below.

        <table border="1">

        <tr><td>%Keyword<td>Type<td>Default<td>Description

        <tr><td><tt>include_q</tt><td>boolean<td>false<td>Some of the
        atoms can be specified as <tt>Q</tt> and given a customizable
        charge.  Such atoms are a point charge that do not have basis
        functions.  If this option is true, then the <tt>Q</tt> atoms are
        included when computing the nuclear charge and the electric field
        due to the nuclear charge.

        <tr><td><tt>include_qq</tt><td>boolean<td>false<td>Some of the
        atoms can be specified as <tt>Q</tt> and given a customizable
        charge.  Such atoms are a point charge that do not have basis
        functions.  If this option is true, then the <tt>Q</tt> atoms are
        included when computing the nuclear repulsion energy and its
        derivatives.

        <tr><td><tt>atominfo</tt><td>AtomInfo<td>library values<td>This
        gives information about each atom, such as the symbol, name, and
        various atomic radii.

        <tr><td><tt>symmetry</tt><td>string<td><tt>C1</tt><td>The
        Schoenflies symbol of the point group.  This is case insensitive.
        It should be a subgroup of D<sub>2h</sub>.  If it is <tt>auto</tt>,
        then the appropriate subgroup of D<sub>2h</sub> will be found.

        <tr><td><tt>symmetry_tolerance</tt><td>double<td>1.0e-4<td>When
        a molecule has symmetry, some atoms may be related by symmetry
        operations.  The distance between given atoms and atoms generated
        by symmetry operations is compared to this threshold to determine
        if they are the same.  If they are the same, then the coordinates
        are cleaned up to make them exactly symmetry equivalent.  If the
        given molecule was produced by a optimization that started in C1
        symmetry, but produced a roughly symmetric structure and you would
        like to begin using symmetry, then this may need to be increased a
        bit to properly symmetrize the molecule.

        <tr><td><tt>symmetry_frame</tt><td>double[3][3]<td>[[1 0 0][0 1
        0][0 0 1]]<td>The symmetry frame.  Ignored for <tt>symmetry =
        auto</tt>.

        <tr><td><tt>origin</tt><td>double[3]<td>[0 0 0]<td>The origin of
        the symmetry frame.  Ignored for <tt>symmetry = auto</tt>.

        <tr><td><tt>redundant_atoms</tt><td>boolean<td>false<td>If true,
        do not generate symmetry equivalent atoms; they are already given
        in the input.  It should not be necessary to specify this option,
        since, by default, if a symmetry operation duplicates an atom, the
        generated atom will not be added to the list of atoms.  Ignored for
        <tt>symmetry = auto</tt>.

        <tr><td><tt>pdb_file</tt><td>string<td>undefined<td>This gives
        the name of a PDB file, from which the nuclear coordinates will be
        read.  If this is given, the following options will be ignored.

        <tr><td><tt>unit</tt><td>string<td>bohr<td>This gives the name
        of the units used for the geometry.  See the Units class for
        information about the known units.  This replaces deprecated
        keywords that are still recognized: <tt>angstrom</tt> and
        <tt>angstroms</tt>.  This is ignored if <tt>pdb_file</tt> is given.

        <tr><td><tt>geometry</tt><td>double[][3]<td>none<td>This gives
        the Cartesian coordinates of the molecule.  This is ignored if
        <tt>pdb_file</tt> is given.

        <tr><td><tt>atoms</tt><td>string[]<td>none<td>This gives the
        chemical element symbol for each atom.  This is ignored if
        <tt>pdb_file</tt> is given.

        <tr><td><tt>ghost</tt><td>boolean[]<td>none<td>If true, the atom
        will be given zero charge.  It will still have basis functions,
        however.  This is used to estimate basis set superposition error.
        This is ignored if <tt>pdb_file</tt> is given.

        <tr><td><tt>charge</tt><td>double[]<td>Z for each atom<td>Allows
        specification of the charge for each atom.  This is ignored if
        <tt>pdb_file</tt> is given.

        <tr><td><tt>atom_labels</tt><td>string[]<td>none<td>This gives a
        user defined atom label for each atom.  This is ignored if
        <tt>pdb_file</tt> is given.

        <tr><td><tt>mass</tt><td>double[]<td>Taken from AtomInfo given by
        the <tt>atominfo</tt> keyword. <td>This gives a user defined mass
        for each atom.  This is ignored if <tt>pdb_file</tt> is given.

        <tr><td><tt>fragment</tt><td>integer[]<td>none<td>Allows to specify
        fragments of Molecule. Fragment indices can be arbitrary
        integers and they do not be consecutive (i.e. one could specify only fragments 1 and 7).
        By default, all atoms belong to fragment 0.
        This feature is relevant only for some computations.
        This keyword is ignored if <tt>pdb_file</tt> is given and its
        values are provided by field resSeq (residue sequence)
        of ATOM or HETATM records.

        </table>

    */
    Molecule(const Ref<KeyVal>&input);

    virtual ~Molecule();

    Molecule& operator=(const Molecule&);

    /// Add an AtomicCenter to the Molecule.
    void add_atom(int Z,double x,double y,double z,
                  const std::string & label = "", double mass = 0.0,
                  int have_charge = 0, double charge = 0.0,
                  int have_fragment = 0, int fragment = 0);

    Atom atom(int i){ return atoms_[i]; }
    const std::vector<Atom>& atoms() const { return atoms_; }

    /// Print information about the molecule.
    virtual void print(std::ostream& =ExEnv::out0()) const;
    virtual void print_parsedkeyval(std::ostream& =ExEnv::out0(),
                                    int print_pg = 1,
                                    int print_unit = 1,
                                    int number_atoms = 1) const;

    /// in which units Molecule was specified and in which units it will be reported
    Ref<Units> geometry_units() const { return geometry_units_; }

    /// Returns the number of atoms in the molecule.
    unsigned int natom() const { return atoms_.size(); }

    int Z(int atom) const { return atoms_[atom].Z(); }
    double &r(int atom, int xyz) { return atoms_[atom].xyz(xyz); }
    const double &r(int atom, int xyz) const { return atoms_[atom].xyz(xyz); }
    const double *r(int atom) const { return atoms_[atom].r(); }
    double mass(int atom) const;
    /** Returns the label explicitly assigned to atom.  If
        no label has been assigned, then null is returned. */
    const char *label(int atom) const;
    /// returns the fragment to which atom belongs to
    int fragment(int atom) const;

    /** Takes an (x, y, z) postion and finds an atom within the
        given tolerance distance.  If no atom is found -1 is returned. */
    int atom_at_position(double *, double tol = 0.05) const;

    /** Returns the index of the atom with the given label.
        If the label cannot be found -1 is returned. */
    int atom_label_to_index(const std::string &label) const;

    /** Returns a vector of the nuclear
        charges of the atoms. */
    std::vector<double> charges() const;

    /// Return the charge of the atom.
    double charge(int iatom) const;

    /// Return true if iatom is a simple point charge
    bool is_Q(int iatom) const;

    /// Returns the total nuclear charge. If include_q is true, this includes
    /// classical charges.
    double total_charge() const;

    /// Returns the sum of atomic numbers of nuclei.
    int total_Z() const;

    /// Sets the PointGroup of the molecule.
    void set_point_group(const Ref<PointGroup>&, double tol=1.0e-7);
    /// Returns the PointGroup of the molecule.
    const Ref<PointGroup>& point_group() const;

    /** Find this molecules true point group (limited to abelian groups).
        If the point group of this molecule is set to the highest point
        group, then the origin must first be set to the center of mass. */
    Ref<PointGroup> highest_point_group(double tol = 1.0e-8) const;

    /** Return 1 if this given axis is a symmetry element for the molecule.
        The direction vector must be a unit vector. */
    int is_axis(SCVector3 &origin,
                SCVector3 &udirection, int order, double tol=1.0e-8) const;

    /** Return 1 if the given plane is a symmetry element for the molecule.
        The perpendicular vector must be a unit vector. */
    int is_plane(SCVector3 &origin, SCVector3 &uperp, double tol=1.0e-8) const;

    /// Return 1 if the molecule has an inversion center.
    int has_inversion(SCVector3 &origin, double tol = 1.0e-8) const;

    /// Returns 1 if the molecule is linear, 0 otherwise.
    int is_linear(double tolerance = 1.0e-5) const;
    /// Returns 1 if the molecule is planar, 0 otherwise.
    int is_planar(double tolerance = 1.0e-5) const;
    /** Sets linear to 1 if the molecular is linear, 0 otherwise.
        Sets planar to 1 if the molecular is planar, 0 otherwise. */
    void is_linear_planar(int&linear,int&planar,double tol = 1.0e-5) const;

    /** Returns a SCVector3 containing the cartesian coordinates of
        the center of mass for the molecule. */
    SCVector3 center_of_mass() const;

    /// Returns the nuclear repulsion energy for the molecule
    double nuclear_repulsion_energy();

    /** Compute the nuclear repulsion energy first derivative with respect
        to the given center. */
    void nuclear_repulsion_1der(int center, double xyz[3]);

    /// Compute the electric field due to the nuclei at the given point.
    void nuclear_efield(const double *position, double* efield);

    /** Compute the electric field due to the given charges at the
        positions of the nuclei at the given point. */
    void nuclear_charge_efield(const double *charges,
                               const double *position, double* efield);

    /** If the molecule contains only symmetry unique atoms, this function
        will generate the other, redundant atoms.  The redundant atom
        will only be generated if there is no other atoms within a distance
        of tol.  If the is another atom and it is not identical, then
        abort will be called. */
    void symmetrize(double tol = 0.5);

    /// Set the point group and then symmetrize.
    void symmetrize(const Ref<PointGroup> &pg, double tol = 0.5);

    /** This will try to carefully correct symmetry errors
        in molecules.  If any atom is out of place by more then
        tol, abort will be called. */
    void cleanup_molecule(double tol = 0.1);

    void translate(const double *r);
    void move_to_com();
    void transform_to_principal_axes(int trans_frame=1);
    void transform_to_symmetry_frame();
    void print_pdb(std::ostream& =ExEnv::out0(), char *title =0) const;

    void read_pdb(const char *filename);

    /** Compute the principal moments of inertia and, possibly, the
        principal axes. */
    void principal_moments_of_inertia(double *evals, double **evecs=0) const;

    /**
     * Return information about symmetry unique and equivalent atoms.
     */
    //@{
    /// Returns the number of symmetry-unique atoms
    int nunique() const { return nuniq_; }
    /// Returns the overall number of the iuniq'th unique atom.
    int unique(int iuniq) const { return equiv_[iuniq][0]; }
    /// Returns the number of atoms equivalent to iuniq.
    int nequivalent(int iuniq) const { return nequiv_[iuniq]; }
    /// Returns the j'th atom equivalent to iuniq.
    int equivalent(int iuniq, int j) const { return equiv_[iuniq][j]; }
    /** Converts an atom number to the number of its generating unique atom.
        The return value is in [0, nunique). */
    int atom_to_unique(int iatom) const { return atom_to_uniq_[iatom]; }
    /** Converts an atom number to the offset of this atom in the list of
        generated atoms. The unique atom itself is allows offset 0.  */
    int atom_to_unique_offset(int iatom) const;
    //@}

    /// Return the number of core electrons.
    int n_core_electrons();

    /// Return the maximum atomic number.
    int max_z();

    /// Return the molecule's AtomInfo object.
    Ref<AtomInfo> atominfo() const { return atominfo_; }

    /// Returns the element name of the atom.
    std::string atom_name(int iatom) const;

    /// Returns the element symbol of the atom.
    std::string atom_symbol(int iatom) const;

    /** If include_q is true, then include the "Q" atoms in the charge and
        efield routines. */
    void set_include_q(bool iq) { include_q_ = iq; }
    /// Returns include_q.  See set_include_q.
    bool include_q() const { return include_q_; }

    /** If include_qq is true, include the coupling between pairs of "Q"
        atoms when computing nuclear repulsion energy and gradients. */
    void set_include_qq(bool iqq) { include_qq_ = iqq; }
    /// Returns include_qq.  See set_include_qq.
    bool include_qq() const { return include_qq_; }

    /// Retrieve the number of "Q" atoms.
    int n_q_atom() const { return q_atoms_.size(); }
    /// Retrieve the "Q" atoms.
    int q_atom(int i) const { return q_atoms_[i]; }

    /// Retrieve the number of non-"Q" atoms.
    int n_non_q_atom() const { return non_q_atoms_.size(); }
    /// Retrieve the of non-"Q" atoms.
    int non_q_atom(int i) const { return non_q_atoms_[i]; }

    /// Determine if any of the atoms have non-standard charge
    bool any_atom_has_charge() const;
    /// Determine if any of the atoms have a fragment label
    bool any_atom_has_fragment() const;
    /// Determine if any of the atoms have a user defined label
    bool any_atom_has_label() const;

    /// returns the origin of the reference coordinate system (the system in which atoms were specified
    /// before the center-of-mass shift
    SCVector3 ref_origin() const { return ref_origin_; }

    void save_data_state(StateOut&);

    virtual ptree& write_xml(ptree& parent, const XMLWriter& writer);

};

/// @}
// end of addtogroup ChemistryMolecule

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
