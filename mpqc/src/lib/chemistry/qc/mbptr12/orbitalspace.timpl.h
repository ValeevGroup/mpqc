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

#ifdef __GNUG__
#pragma interface
#endif

#ifndef _mpqc_src_lib_chemistry_qc_mbptr12_orbitalspacetimpl_h
#define _mpqc_src_lib_chemistry_qc_mbptr12_orbitalspacetimpl_h

// includes go here
#include <chemistry/qc/mbptr12/orbitalspace.h>
#include <algorithm>
#include <cassert>
#include <vector>

namespace sc {

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
    assert(norbs != 0);
    const unsigned int max_orbsym = *(std::max_element(orbsyms.begin(), orbsyms.end()));
    const unsigned int min_orbsym = *(std::min_element(orbsyms.begin(), orbsyms.end()));
    const unsigned int nirreps = basis->molecule()->point_group()->char_table().order();
    if (min_orbsym >= nirreps || max_orbsym >= nirreps)
      throw ProgrammingError("OrderedOrbitalSpace -- invalid orbital symmetry array",__FILE__,__LINE__);

    // construct vector of MolecularOrbital objects
    std::vector<MolecularOrbital> orbs;
    for(size_t o=0; o<norbs; ++o) {
      using detail::MolecularOrbitalAttributes;
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

} // end of namespace sc

#endif // end of header guard

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
