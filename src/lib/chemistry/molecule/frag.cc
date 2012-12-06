//
// frag.cc
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

#include <algorithm>
#include <cassert>
#include <chemistry/molecule/frag.h>
#include <util/misc/scexception.h>
#include <chemistry/molecule/coor.h>
#include <util/state/statein.h>
#include <util/state/stateout.h>
#include <util/container/bitarray.h>

using namespace sc;

ClassDesc
MolecularFragment::class_desc_(typeid(MolecularFragment),
                               "MolecularFragment",
                               1,               // version
                               "public Molecule", // must match parent
                               0,               // change to create<MolecularFragment> if this class is DefaultConstructible
                               create<MolecularFragment>, // change to 0 if this class is not KeyValConstructible
                               create<MolecularFragment>  // change to 0 if this class is not StateInConstructible
);

MolecularFragment::MolecularFragment(const Ref<KeyVal>& kv) :
    Molecule()
{
  std::set<int> atoms_in_frag;
  process_keyval(kv, protomol_, fragment_, includes_atoms_, excludes_atoms_, atoms_in_frag);

  // need to include the rest of atoms as ghosts?
  const bool add_ghosts = kv->booleanvalue("add_ghosts", KeyValValueboolean(false));

  // add atoms to the base molecule
  const int natoms_in_protomol = protomol_->natom();
  for(int atom=0; atom<natoms_in_protomol; ++atom) {
    // add as ghost?
    const bool add_to_fragment = (atoms_in_frag.find(atom) != atoms_in_frag.end());
    if (add_to_fragment == false && add_ghosts == false)
      continue;
    const double Z = protomol_->Z(atom);
    const SCVector3 O = protomol_->ref_origin();
    double r[3]; for(int xyz=0; xyz<3; ++xyz) r[xyz] = protomol_->r(atom, xyz);
    const char* clabel = protomol_->label(atom);
    const std::string label(clabel == 0 ? "" : clabel);
    const double mass = protomol_->mass(atom);
    const int have_charge = 1;
    const double charge = add_to_fragment ? protomol_->charge(atom) : 0.0;
    // keep fragment info
    const int have_fragment = 1;
    const int fragment = protomol_->fragment(atom);
    Molecule::add_atom(Z, r[0], r[1], r[2], label, mass,
                       have_charge, charge,
                       have_fragment, fragment);
  }
  for(int xyz=0; xyz<3; ++xyz) ref_origin_[xyz] = protomol_->ref_origin()[xyz];

  // determine point group
  const std::string symmetry = kv->stringvalue("symmetry", KeyValValuestring("auto"));
  double symtol = kv->doublevalue("symmetry_tolerance",
                                  KeyValValuedouble(1.0e-4));
  Ref<PointGroup> pg;
  if (symmetry == "auto") {
    pg = highest_point_group(symtol);
  }
  else {
    pg = new PointGroup(kv);

    // origin provided by user in units of protomol, convert to atomic units
    const double conv = protomol_->geometry_units()->to_atomic_units();
    for (int i=0; i<3; i++) {
      pg->origin()[i] *= conv;
    }
  }
  set_point_group(pg, symtol*10.0);

  // OCD: report fragment coordinates in same units as protomol
  geometry_units_ = protomol_->geometry_units();

  protomol_->print_parsedkeyval();
  this->print_parsedkeyval();
}

MolecularFragment::MolecularFragment(StateIn& si) : Molecule(si) {
  protomol_ << SavableState::restore_state(si);
  si.get(fragment_);
  si.get(includes_atoms_);
  si.get(excludes_atoms_);
}

MolecularFragment::~MolecularFragment() {
  // this may be necessary if this is a templated class
  const bool make_sure_class_desc_initialized = (&class_desc_ != 0);
}

void
MolecularFragment::save_data_state(StateOut& so) {
  Molecule::save_data_state(so);
  SavableState::save_state(protomol_.pointer(), so);
  so.put(fragment_);
  so.put(includes_atoms_);
  so.put(excludes_atoms_);
}

void
MolecularFragment::process_keyval(const Ref<KeyVal>& kv,
                                  Ref<Molecule>& protomol,
                                  int& fragment,
                                  std::set<int>& includes_atoms,
                                  std::set<int>& excludes_atoms,
                                  std::set<int>& atoms_in_frag)
{
  typedef std::set<int>::const_iterator citer;

  protomol << kv->describedclassvalue("molecule");
  if (protomol.null()) throw InputError("missing keyword", __FILE__, __LINE__,
                                        "MolecularFragment::molecule", "");
  const int natoms = protomol->natom();

  // filter the atoms based on fragment
  fragment = kv->intvalue("fragment", KeyValValueint(0));
  {
    for(int atom=0; atom<natoms; ++atom)
      if (protomol->fragment(atom) == fragment)
        atoms_in_frag.insert(atom);
  }

  // compute the adjacency matrix and the subgraphs of the molecular graph
  const BitArrayLTri adjmat = IntCoorGen::adjacency_matrix(*protomol, 1.1);
  std::vector<std::set<int> > subgraphs = IntCoorGen::find_disconnected_subgraphs(adjmat);
  const size_t nsubgraphs = subgraphs.size();

  // filter the atoms based on connectivity: include only the atoms connected to the ones given by includes_atoms
  // skip this all atoms are connected
  if (nsubgraphs > 1) {
    if (kv->exists("includes_atoms")) {
      Keyword kw(kv, "includes_atoms");
      kw >> includes_atoms;

      std::set<int> subgraphs_to_include;
      // for each atom in the includes_atoms list that is found in atoms_in_frag record the subgraph to which it belongs
      for(citer i=includes_atoms.begin(); i!=includes_atoms.end(); ++i) {
        if (atoms_in_frag.find(*i) != atoms_in_frag.end()) {
          // which subgraph?
          for(unsigned int s=0; s<nsubgraphs; ++s) {
            if (subgraphs[s].find(*i) != subgraphs[s].end()) {
              subgraphs_to_include.insert(s);
              break;
            }
          }
        }
      }

      // make union of the selected subgraphs
      std::set<int> subgraph_union;
      for(int s=0; s<subgraphs_to_include.size(); ++s) { // guaranteed to have more than subgraph
        std::set<int> new_subgraph_union;
        std::set_union(subgraph_union.begin(), subgraph_union.end(),
                       subgraphs[s].begin(), subgraphs[s].end(),
                       inserter(new_subgraph_union, new_subgraph_union.end())
                       );
        std::swap(subgraph_union,new_subgraph_union);
      }

      // intersection of subgraph_union and atoms_in_frag gives new atoms_in_frag
      std::set<int> new_atoms_in_frag;
      std::set_intersection(atoms_in_frag.begin(), atoms_in_frag.end(),
                            subgraph_union.begin(), subgraph_union.end(),
                            inserter(new_atoms_in_frag, new_atoms_in_frag.begin())
                           );
      std::swap(atoms_in_frag, new_atoms_in_frag);
    }
  }

  // filter the atoms based on connectivity: include only the atoms NOT connected to the ones given by excludes_atoms
  if (kv->exists("excludes_atoms")) {
    Keyword kw(kv, "excludes_atoms");
    kw >> excludes_atoms;

    std::set<int> subgraphs_to_exclude;
    // for each atom in the excludes_atoms list that is found in atoms_in_frag record the subgraph to which it belongs
    for(citer i=excludes_atoms.begin(); i!=excludes_atoms.end(); ++i) {
      if (atoms_in_frag.find(*i) != atoms_in_frag.end()) {
        // which subgraph?
        for(unsigned int s=0; s<nsubgraphs; ++s) {
          if (subgraphs[s].find(*i) != subgraphs[s].end()) {
            subgraphs_to_exclude.insert(s);
            break;
          }
        }
      }
    }

    // make union of the selected subgraphs
    std::set<int> subgraph_union;
    for(int s=0; s<subgraphs_to_exclude.size(); ++s) { // guaranteed to have more than subgraph
      std::set<int> new_subgraph_union;
      std::set_union(subgraph_union.begin(), subgraph_union.end(),
                     subgraphs[s].begin(), subgraphs[s].end(),
                     inserter(new_subgraph_union, new_subgraph_union.end())
                     );
      std::swap(subgraph_union,new_subgraph_union);
    }

    // set-theoretic difference atoms_in_frag - subgraph_union gives new atoms_in_frag
    std::set<int> new_atoms_in_frag;
    std::set_difference(atoms_in_frag.begin(), atoms_in_frag.end(),
                        subgraph_union.begin(), subgraph_union.end(),
                        inserter(new_atoms_in_frag, new_atoms_in_frag.begin())
                       );
    std::swap(atoms_in_frag, new_atoms_in_frag);
  }

  // have any atoms left?
  if (atoms_in_frag.empty())
    throw InputError("no atoms in fragment", __FILE__, __LINE__);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
