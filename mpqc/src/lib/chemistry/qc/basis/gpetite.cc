//
// gpetite.cc --- implementation of GPetite4 and helpers
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
// Maintainer: SNL
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <util/misc/formio.h>

#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/gpetite.h>

using namespace std;
using namespace sc;

////////////////////////////////////////////////////////////////////////////

canonical_aaaa::canonical_aaaa()
{
}

canonical_aaaa::canonical_aaaa(const Ref<GaussianBasisSet> bi,
                               const Ref<GaussianBasisSet> bj,
                               const Ref<GaussianBasisSet> bk,
                               const Ref<GaussianBasisSet> bl)
{
}

////////////////////////////////////////////////////////////////////////////

canonical_aabc::canonical_aabc(const Ref<GaussianBasisSet> bi,
                               const Ref<GaussianBasisSet> bj,
                               const Ref<GaussianBasisSet> bk,
                               const Ref<GaussianBasisSet> bl)
{
  nk_ = bk->nshell();
  nl_ = bl->nshell();
}

////////////////////////////////////////////////////////////////////////////


canonical_abcc::canonical_abcc(const Ref<GaussianBasisSet> bi,
                               const Ref<GaussianBasisSet> bj,
                               const Ref<GaussianBasisSet> bk,
                               const Ref<GaussianBasisSet> bl)
{
  nj_ = bj->nshell();
  int nk = bk->nshell();
  nkl_ = (nk*long(nk+1))>>1;
}

////////////////////////////////////////////////////////////////////////////


canonical_aabb::canonical_aabb(const Ref<GaussianBasisSet> bi,
                               const Ref<GaussianBasisSet> bj,
                               const Ref<GaussianBasisSet> bk,
                               const Ref<GaussianBasisSet> bl)
{
  int ni = bi->nshell();
  nij_ = (ni*long(ni+1))>>1;
}

////////////////////////////////////////////////////////////////////////////


canonical_abab::canonical_abab(const Ref<GaussianBasisSet> bi,
                               const Ref<GaussianBasisSet> bj,
                               const Ref<GaussianBasisSet> bk,
                               const Ref<GaussianBasisSet> bl)
{
  nj_ = bj->nshell();
}

////////////////////////////////////////////////////////////////////////////

canonical_abcd::canonical_abcd(const Ref<GaussianBasisSet> bi,
                               const Ref<GaussianBasisSet> bj,
                               const Ref<GaussianBasisSet> bk,
                               const Ref<GaussianBasisSet> bl)
{
  ni_ = bi->nshell();
  nj_ = bj->nshell();
  nk_ = bk->nshell();
}


/////////////////////////////////////////////////////////////////////////////

GenPetite4::GenPetite4(const Ref<GaussianBasisSet> &b1,
                       const Ref<GaussianBasisSet> &b2,
                       const Ref<GaussianBasisSet> &b3,
                       const Ref<GaussianBasisSet> &b4)
{
  int **atom_map;
  b1_ = b1;
  b2_ = b2;
  b3_ = b3;
  b4_ = b4;

  ng_ = b1->molecule()->point_group()->char_table().order();
  if (b2->molecule()->point_group()->char_table().order() != ng_
      || b3->molecule()->point_group()->char_table().order() != ng_
      || b4->molecule()->point_group()->char_table().order() != ng_) {
    throw std::runtime_error("GPetite4: not all point groups are the same");
  }
  c1_ =  (ng_ == 1);

  atom_map = compute_atom_map(b1);
  shell_map_i_ = compute_shell_map(atom_map,b1);
  delete_atom_map(atom_map,b1);

  atom_map = compute_atom_map(b2);
  shell_map_j_ = compute_shell_map(atom_map,b2);
  delete_atom_map(atom_map,b2);

  atom_map = compute_atom_map(b3);
  shell_map_k_ = compute_shell_map(atom_map,b3);
  delete_atom_map(atom_map,b3);

  atom_map = compute_atom_map(b4);
  shell_map_l_ = compute_shell_map(atom_map,b4);
  delete_atom_map(atom_map,b4);
}

GenPetite4::~GenPetite4() {
  delete_shell_map(shell_map_i_,b1_);
  delete_shell_map(shell_map_j_,b2_);
  delete_shell_map(shell_map_k_,b3_);
  delete_shell_map(shell_map_l_,b4_);
}

/////////////////////////////////////////////////////////////////////////////

template <class C4> GPetite4<C4>::GPetite4(const Ref<GaussianBasisSet> &b1,
         const Ref<GaussianBasisSet> &b2,
         const Ref<GaussianBasisSet> &b3,
         const Ref<GaussianBasisSet> &b4,
         const C4& c): GenPetite4(b1,b2,b3,b4), c_(c)
{
}

template <class C4> GPetite4<C4>::~GPetite4()
{
}

/////////////////////////////////////////////////////////////////////////////

Ref<GenPetite4>
sc::construct_gpetite(const Ref<GaussianBasisSet> &b1,
                  const Ref<GaussianBasisSet> &b2,
                  const Ref<GaussianBasisSet> &b3,
                  const Ref<GaussianBasisSet> &b4)
{
  if (b1 == b2 && b1 == b3 && b1 == b4) {
    canonical_aaaa c4(b1,b2,b3,b4);
    return new GPetite4<canonical_aaaa>(b1,b2,b3,b4,c4);
  }
  else if (b1 == b2 && b3 != b4) {
    canonical_aabc c4(b1,b2,b3,b4);
    return new GPetite4<canonical_aabc>(b1,b2,b3,b4,c4);
  }
  else if (b1 != b2 && b3 == b4) {
    canonical_abcc c4(b1,b2,b3,b4);
    return new GPetite4<canonical_abcc>(b1,b2,b3,b4,c4);
  }
  else if (b1 == b3 && b2 == b4) {
    canonical_abab c4(b1,b2,b3,b4);
    return new GPetite4<canonical_abab>(b1,b2,b3,b4,c4);
  }
  else if (b1 == b2 && b3 == b4) {
    canonical_aabb c4(b1,b2,b3,b4);
    return new GPetite4<canonical_aabb>(b1,b2,b3,b4,c4);
  }
  else {
    canonical_abcd c4(b1,b2,b3,b4);
    return new GPetite4<canonical_abcd>(b1,b2,b3,b4,c4);
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
