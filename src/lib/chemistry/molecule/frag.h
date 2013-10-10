//
// frag.h
//
// Copyright (C) 2011 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
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

#ifndef _mpqc_src_lib_chemistry_molecule_frag_h
#define _mpqc_src_lib_chemistry_molecule_frag_h

#include <chemistry/molecule/molecule.h>

namespace sc {

  /// @addtogroup ChemistryMolecule
  /// @{

  /// MolecularFragment is a Molecule that is a fragment of another Molecule object
  class MolecularFragment : public Molecule {
    public:
      /** A KeyVal constructor is used to generate a MolecularFragment
          object from the input. The full list of keywords
          that are accepted is below.

          <table border="1">

          <tr><td><b>%Keyword</b><td><b>Type</b><td><b>Default</b><td><b>Description</b>

          <tr><td><tt>molecule</tt><td>Molecule<td>none<td>the Molecule object that this molecule is part of (proto-molecule).

          <tr><td><tt>fragment</tt><td>integer<td>0<td>Specifies atoms of which fragment defined in the proto-molecule
          object will be considered for including into this molecule (see Molecule KeyVal documentation). By default,
          all atoms will be considered.

          <tr><td><tt>includes_atoms</tt><td>integer[] <td>none<td>(optional) only atoms connected by
          "normal" chemical bonds to the atoms specified by this keyword will be considered (atomic indices are 0-based).
          The normal chemical bonds are generated for the proto-molecule object using
          the same distance criterion as that used by IntCoorGen object
          to generate bonds (see IntCoorGen documentation).

          <tr><td><tt>excludes_atoms</tt><td>integer[] <td>none<td>(optional) only atoms NOT connected by
          "normal" chemical bonds to the atoms specified by this keyword will be considered. See documentation for <tt>includes_atoms</tt>.

        <tr><td><tt>symmetry</tt><td>string<td><tt>auto</tt><td>The
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
        the symmetry frame. Specified in the same units as used by proto-molecule.
        Ignored for <tt>symmetry = auto</tt>.

        <tr><td><tt>add_ghosts</tt><td>boolean<td>false<td>If true, the atoms
        in proto-molecule that are not of this fragment will be included as ghosts
        (zero-charge atoms that carry basis functions; see documentation for
        <tt>ghost</tt> keyword for Molecule).

          </table>
          Note that although MolecularFragment is derived from Molecule, its KeyVal constructor does not
          accept all Molecule keywords.
       */
      MolecularFragment(const Ref<KeyVal>& kv);
      MolecularFragment(StateIn&);
      virtual ~MolecularFragment();
      void save_data_state(StateOut&);

    private:
      static ClassDesc class_desc_;
      Ref<Molecule> protomol_;
      int fragment_;
      std::set<int> includes_atoms_;
      std::set<int> excludes_atoms_;

      static void process_keyval(const Ref<KeyVal>& kv,
                                 Ref<Molecule>& protomol,
                                 int& fragment,
                                 std::set<int>& includes_atoms,
                                 std::set<int>& excludes_atoms,
                                 std::set<int>& atoms_in_frag);
  };

  /// @}
  // end of addtogroup ChemistryMolecule

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
