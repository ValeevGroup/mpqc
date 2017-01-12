//
// orbitalspace.timpl.h
//
// Copyright (C) 2009 Edward Valeev
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

#ifndef _mpqc_src_lib_chemistry_qc_mbptr12_orbitalspacetimpl_h
#define _mpqc_src_lib_chemistry_qc_mbptr12_orbitalspacetimpl_h

// includes go here
#include <chemistry/qc/wfn/orbitalspace.h>
#include <algorithm>
#include <cassert>
#include <vector>

namespace sc {

  template <typename Order>
  ClassDesc
  OrderedOrbitalSpace<Order>::class_desc_(typeid(this_type),
                                          (std::string("OrderedOrbitalSpace<") +
                                           std::string(typeid(Order).name()) +
                                           std::string(">")
                                          ).c_str(), 1,
                                          "public OrbitalSpace", 0, 0,
                                          create<this_type> );

  template <typename Order>
  OrderedOrbitalSpace<Order>::OrderedOrbitalSpace(StateIn& si) :
    OrbitalSpace(si) {}

  template <typename Order>
  void
  OrderedOrbitalSpace<Order>::save_data_state(StateOut& so) {
    OrbitalSpace::save_data_state(so);
  }

  template <typename Order>
  OrderedOrbitalSpace<Order>::~OrderedOrbitalSpace() {
  }

  template <typename Order>
  OrderedOrbitalSpace<Order>::OrderedOrbitalSpace(const std::string& id,
                                           const std::string& name,
                                           const Ref<GaussianBasisSet>& basis,
                                           const Ref<Integral>& integral,
                                           const RefSCMatrix& coefs,
                                           const RefDiagSCMatrix& evals,
                                           const RefDiagSCMatrix& occnums,
                                           const std::vector<unsigned int>& orbsyms,
                                           const Order& order) :
    OrbitalSpace() {

    // validate input
    const size_t norbs = coefs.coldim().n();
    if (orbsyms.empty() == false) {
      const unsigned int max_orbsym = *(std::max_element(orbsyms.begin(), orbsyms.end()));
      const unsigned int min_orbsym = *(std::min_element(orbsyms.begin(), orbsyms.end()));
      const unsigned int nirreps = basis->molecule()->point_group()->char_table().order();
      if (min_orbsym >= nirreps || max_orbsym >= nirreps)
        throw ProgrammingError("OrderedOrbitalSpace -- invalid orbital symmetry array",__FILE__,__LINE__);
    }

    // construct vector of MolecularOrbital objects
    std::vector<MolecularOrbital> orbs;
    for(size_t o=0; o<norbs; ++o) {
      orbs.push_back(MolecularOrbital(o,
                                      MolecularOrbitalAttributes(orbsyms.at(o),
                                                                 evals.get_element(o),
                                                                 occnums.get_element(o)
                                      )
      ));
    }

    // sort
    std::stable_sort(orbs.begin(), orbs.end(), order);

    // convert vector of MolecularOrbitals to BlockedOrbitals
    std::vector<BlockedOrbital> blocked_orbs;
    for(size_t o=0; o<norbs; ++o) {
      blocked_orbs.push_back(BlockedOrbital(orbs[o].index(),
                                            order.block(orbs[o])
      ));
    }

    init(id, name, basis, integral, coefs, evals, orbsyms, order.nblocks(), blocked_orbs);

  }

  //////////////////////////////////////////////////////////

  template <typename Order>
  ClassDesc
  OrderedSpinOrbitalSpace<Order>::class_desc_(typeid(this_type),
                                          (std::string("OrderedSpinOrbitalSpace<") +
                                           std::string(typeid(Order).name()) +
                                           std::string(">")
                                          ).c_str(), 1,
                                          "public OrbitalSpace", 0, 0,
                                          create<this_type> );

  template <typename Order>
  OrderedSpinOrbitalSpace<Order>::OrderedSpinOrbitalSpace(StateIn& si) :
    OrbitalSpace(si) {}

  template <typename Order>
  void
  OrderedSpinOrbitalSpace<Order>::save_data_state(StateOut& so) {
    OrbitalSpace::save_data_state(so);
  }

  template <typename Order>
  OrderedSpinOrbitalSpace<Order>::~OrderedSpinOrbitalSpace() {
  }

  template <typename Order>
  OrderedSpinOrbitalSpace<Order>::OrderedSpinOrbitalSpace(const std::string& id,
                                           const std::string& name,
                                           const Ref<GaussianBasisSet>& basis,
                                           const Ref<Integral>& integral,
                                           const RefSCMatrix& coefs_a,
                                           const RefSCMatrix& coefs_b,
                                           const RefDiagSCMatrix& evals_a,
                                           const RefDiagSCMatrix& evals_b,
                                           const RefDiagSCMatrix& occnums_a,
                                           const RefDiagSCMatrix& occnums_b,
                                           const std::vector<unsigned int>& orbsyms_a,
                                           const std::vector<unsigned int>& orbsyms_b,
                                           const Order& order) :
    OrbitalSpace() {

    // validate input
    const size_t norbs = coefs_a.coldim().n();
    MPQC_ASSERT(norbs != 0);
    MPQC_ASSERT(norbs == coefs_b.coldim().n());
    const unsigned int max_orbsym_a = *(std::max_element(orbsyms_a.begin(), orbsyms_a.end()));
    const unsigned int min_orbsym_a = *(std::min_element(orbsyms_a.begin(), orbsyms_a.end()));
    const unsigned int max_orbsym_b = *(std::max_element(orbsyms_b.begin(), orbsyms_b.end()));
    const unsigned int min_orbsym_b = *(std::min_element(orbsyms_b.begin(), orbsyms_b.end()));
    const unsigned int max_orbsym = std::max(max_orbsym_a,max_orbsym_b);
    const unsigned int min_orbsym = std::min(min_orbsym_a,min_orbsym_b);
    const unsigned int nirreps = basis->molecule()->point_group()->char_table().order();
    if (min_orbsym >= nirreps || max_orbsym >= nirreps)
      throw ProgrammingError("OrderedSpinOrbitalSpace -- invalid orbital symmetry arrays",__FILE__,__LINE__);

    /////////////
    // Merge alpha and beta orbitals:
    /////////////

    // 1) construct vector of MolecularOrbital objects
    std::vector<MolecularSpinOrbital> orbs;
    for(size_t o=0; o<norbs; ++o) { // alpha spin
      orbs.push_back(MolecularSpinOrbital(o,
                                      MolecularSpinOrbitalAttributes(orbsyms_a.at(o),
                                                                     evals_a.get_element(o),
                                                                     occnums_a.get_element(o),
                                                                     Alpha
                                      )
      ));
    }
    for(size_t o=0; o<norbs; ++o) { // beta spin
      orbs.push_back(MolecularSpinOrbital(o + norbs,
                                      MolecularSpinOrbitalAttributes(orbsyms_b.at(o),
                                                                     evals_b.get_element(o),
                                                                     occnums_b.get_element(o),
                                                                     Beta
                                      )
      ));
    }

    // sort
    std::stable_sort(orbs.begin(), orbs.end(), order);

    // convert vector of MolecularOrbitals to BlockedOrbitals
    std::vector<BlockedOrbital> blocked_orbs;
    for(size_t o=0; o<norbs*2; ++o) {
      const unsigned index = orbs[o].index();
      const unsigned block = order.block(orbs[o]);
      BlockedOrbital orb(index,block);
      blocked_orbs.push_back(orb);
    }

    // 2) merge coefficients, eigenvalues, and occupations
    RefSCDimension orbdim;
    {
      const unsigned int nblocks = coefs_a.coldim()->blocks()->nblock() * 2;
      // build new blocked dimension
      int* nfunc_per_block = new int[nblocks];
      for (unsigned int i = 0; i < nblocks/2; ++i)
        nfunc_per_block[i] = coefs_a.coldim()->blocks()->size(i);
      for (unsigned int i = 0, ii=nblocks/2; i < nblocks/2; ++i, ++ii)
        nfunc_per_block[ii] = coefs_a.coldim()->blocks()->size(i);
      orbdim = new SCDimension(norbs * 2, nblocks, nfunc_per_block, id.c_str());
      if (norbs) {
        for (unsigned int i = 0; i < nblocks; ++i)
          orbdim->blocks()->set_subdim(i, new SCDimension(nfunc_per_block[i]));
      }
      delete[] nfunc_per_block;
    }
    RefSCMatrix coefs = coefs_a.kit()->matrix(coefs_a.rowdim(),orbdim);
    RefDiagSCMatrix evals = evals_a.kit()->diagmatrix(orbdim);
    std::vector<unsigned int> orbsyms(norbs*2);
    const unsigned int nao = coefs_a.rowdim().n();
    for (unsigned int i = 0, ii=0; i < norbs; ++i, ++ii) { // alpha
      for(unsigned int ao=0; ao<nao; ++ao) {
        coefs(ao,ii) = coefs_a(ao,i);
      }
      evals(ii) = evals_a(i);
      orbsyms[ii] = orbsyms_a[i];
    }
    for (unsigned int i = 0, ii=norbs; i < norbs; ++i, ++ii) { // alpha
      for(unsigned int ao=0; ao<nao; ++ao) {
        coefs(ao,ii) = coefs_b(ao,i);
      }
      evals(ii) = evals_b(i);
      orbsyms[ii] = orbsyms_b[i];
    }

    init(id, name, basis, integral, coefs, evals, orbsyms, order.nblocks(), blocked_orbs);

  }


} // end of namespace sc

#endif // end of header guard

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
