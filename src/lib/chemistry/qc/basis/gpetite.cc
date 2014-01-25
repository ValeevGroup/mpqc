//
// gpetite.cc --- implementation of GenericPetiteList4 and helpers
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

#include <util/misc/formio.h>
#include <util/misc/scexception.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/gpetite.h>

using namespace std;
using namespace sc;

////////////////////////////////////////////////////////////////////////////

canonical_aaaa::canonical_aaaa()
{
}

canonical_aaaa::canonical_aaaa(const Ref<GaussianBasisSet>& bi,
                               const Ref<GaussianBasisSet>& bj,
                               const Ref<GaussianBasisSet>& bk,
                               const Ref<GaussianBasisSet>& bl)
{
}

////////////////////////////////////////////////////////////////////////////

canonical_aabc::canonical_aabc(const Ref<GaussianBasisSet>& bi,
                               const Ref<GaussianBasisSet>& bj,
                               const Ref<GaussianBasisSet>& bk,
                               const Ref<GaussianBasisSet>& bl)
{
  nk_ = bk->nshell();
  nl_ = bl->nshell();
}

////////////////////////////////////////////////////////////////////////////


canonical_abcc::canonical_abcc(const Ref<GaussianBasisSet>& bi,
                               const Ref<GaussianBasisSet>& bj,
                               const Ref<GaussianBasisSet>& bk,
                               const Ref<GaussianBasisSet>& bl)
{
  nj_ = bj->nshell();
  int nk = bk->nshell();
  nkl_ = (nk*long(nk+1))>>1;
}

////////////////////////////////////////////////////////////////////////////


canonical_aabb::canonical_aabb(const Ref<GaussianBasisSet>& bi,
                               const Ref<GaussianBasisSet>& bj,
                               const Ref<GaussianBasisSet>& bk,
                               const Ref<GaussianBasisSet>& bl)
{
  int ni = bi->nshell();
  nij_ = (ni*long(ni+1))>>1;
}

////////////////////////////////////////////////////////////////////////////


canonical_abab::canonical_abab(const Ref<GaussianBasisSet>& bi,
                               const Ref<GaussianBasisSet>& bj,
                               const Ref<GaussianBasisSet>& bk,
                               const Ref<GaussianBasisSet>& bl)
{
  nj_ = bj->nshell();
}

////////////////////////////////////////////////////////////////////////////

canonical_abcd::canonical_abcd(const Ref<GaussianBasisSet>& bi,
                               const Ref<GaussianBasisSet>& bj,
                               const Ref<GaussianBasisSet>& bk,
                               const Ref<GaussianBasisSet>& bl)
{
  ni_ = bi->nshell();
  nj_ = bj->nshell();
  nk_ = bk->nshell();
}


////////////////////////////////////////////////////////////////////////////

canonical_aa::canonical_aa(const Ref<GaussianBasisSet>& bi,
                           const Ref<GaussianBasisSet>& bj)
{
}

////////////////////////////////////////////////////////////////////////////


canonical_ab::canonical_ab(const Ref<GaussianBasisSet>& bi,
                           const Ref<GaussianBasisSet>& bj)
{
  ni_ = bi->nshell();
}

/////////////////////////////////////////////////////////////////////////////

GPetiteList4::GPetiteList4(const Ref<GaussianBasisSet> &b1,
                       const Ref<GaussianBasisSet> &b2,
                       const Ref<GaussianBasisSet> &b3,
                       const Ref<GaussianBasisSet> &b4)
{
  int **atom_map;
  b1_ = b1;
  b2_ = b2;
  b3_ = b3;
  b4_ = b4;

  const Ref<PointGroup>& pg = b1->molecule()->point_group();
  ng_ = pg->char_table().order();
  if (not b2->molecule()->point_group()->equiv(pg)
      || not b3->molecule()->point_group()->equiv(pg)
      || not b4->molecule()->point_group()->equiv(pg)) {
    throw ProgrammingError("GenericPetiteList4: not all point groups are the same",__FILE__,__LINE__);
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

GPetiteList4::~GPetiteList4() {
  delete_shell_map(shell_map_i_,b1_);
  delete_shell_map(shell_map_j_,b2_);
  delete_shell_map(shell_map_k_,b3_);
  delete_shell_map(shell_map_l_,b4_);
}

/////////////////////////////////////////////////////////////////////////////

template <class C4> GenericPetiteList4<C4>::GenericPetiteList4(const Ref<GaussianBasisSet> &b1,
         const Ref<GaussianBasisSet> &b2,
         const Ref<GaussianBasisSet> &b3,
         const Ref<GaussianBasisSet> &b4):
           GPetiteList4(b1,b2,b3,b4),
           c_(b1,b2,b3,b4)
{
}

template <class C4> GenericPetiteList4<C4>::~GenericPetiteList4()
{
}

/////////////////////////////////////////////////////////////////////////////

GPetiteList2::GPetiteList2(const Ref<GaussianBasisSet> &b1,
                           const Ref<GaussianBasisSet> &b2)
{
  int **atom_map;
  b1_ = b1;
  b2_ = b2;

  ng_ = b1->molecule()->point_group()->char_table().order();
  if (b2->molecule()->point_group()->char_table().order() != ng_) {
    throw ProgrammingError("GenericPetiteList2: not all point groups are the same",__FILE__,__LINE__);
  }
  c1_ =  (ng_ == 1);

  atom_map = compute_atom_map(b1);
  shell_map_i_ = compute_shell_map(atom_map,b1);
  delete_atom_map(atom_map,b1);

  atom_map = compute_atom_map(b2);
  shell_map_j_ = compute_shell_map(atom_map,b2);
  delete_atom_map(atom_map,b2);
}

GPetiteList2::~GPetiteList2() {
  delete_shell_map(shell_map_i_,b1_);
  delete_shell_map(shell_map_j_,b2_);
}

/////////////////////////////////////////////////////////////////////////////

template <class C2> GenericPetiteList2<C2>::GenericPetiteList2(const Ref<GaussianBasisSet> &b1,
                                                               const Ref<GaussianBasisSet> &b2):
           GPetiteList2(b1,b2),
           c_(b1,b2)
{
}

template <class C2> GenericPetiteList2<C2>::~GenericPetiteList2()
{
}

void
sc::symmetrize(const Ref<GPetiteList2>& plist12,
               const Ref<Integral>& integral,
               const RefSymmSCMatrix& skel,
               const RefSymmSCMatrix& sym)
{
  // validate input
  if (plist12->basis1()->equiv(plist12->basis2()) != true)
    throw ProgrammingError("sc::symmetrize(RefSymmSCMatrix) -- petite list's basis sets are not equivalent",__FILE__,__LINE__);
  if (skel.dim().n() != plist12->basis1()->nbasis())
    throw ProgrammingError("sc::symmetrize(RefSymmSCMatrix) -- dimension of skel doesn't match the number of basis functions",__FILE__,__LINE__);
  if (skel.dim().n() != sym.dim().n())
    throw ProgrammingError("sc::symmetrize(RefSymmSCMatrix) -- dimensions of skel and sym don't match",__FILE__,__LINE__);

  const Ref<GaussianBasisSet>& bs1 = plist12->basis1();
  integral->set_basis(bs1);
  Ref<PetiteList> pl1 = integral->petite_list();

  // SO basis is always blocked, so first make sure skel is blocked
  RefSymmSCMatrix bskel = dynamic_cast<BlockedSymmSCMatrix*>(skel.pointer());
  if (bskel.null()) {
    bskel = bs1->so_matrixkit()->symmmatrix(pl1->AO_basisdim());
    bskel->convert(skel);
  }

  // if C1, then do nothing
  CharacterTable ct = plist12->point_group()->char_table();
  if (ct.order() == 1) {
    sym.assign(bskel);
    return;
  }

  RefSCMatrix aoso = pl1->aotoso();
  BlockedSCMatrix *lu = dynamic_cast<BlockedSCMatrix*>(aoso.pointer());
  for (int b=0; b < lu->nblocks(); b++) {
    if (lu->block(b).null())
      continue;

    const int ir = ct.which_irrep(b);

    double skal = (double)ct.order()/(double)ct.gamma(ir).degeneracy();
    skal = sqrt(skal);
    lu->block(b).scale(skal);
  }

  sym.assign(0.0);
  sym.accumulate_transform(aoso,bskel,SCMatrix::TransposeTransform);
  aoso=0;

  // loop through blocks and finish symmetrizing degenerate blocks
  BlockedSymmSCMatrix *la = dynamic_cast<BlockedSymmSCMatrix*>(sym.pointer());
  for (int b=0; b < la->nblocks(); b++) {
    if (la->block(b).null())
      continue;

    const int ir=ct.which_irrep(b);

    if (ct.gamma(ir).degeneracy()==1)
      continue;

    if (ct.gamma(ir).complex()) {
      const int nbf = pl1->nfunction(ir)/2;

      RefSymmSCMatrix irrep = la->block(b).get_subblock(0, nbf-1);
      irrep.accumulate(la->block(b).get_subblock(nbf, 2*nbf-1));
      la->block(b).assign_subblock(irrep,  0, nbf-1);
      la->block(b).assign_subblock(irrep,nbf, 2*nbf-1);

      RefSCMatrix sub = la->block(b).get_subblock(nbf, 2*nbf-1, 0, nbf-1);
      RefSCMatrix subt = sub.t();
      subt.scale(-1.0);
      sub.accumulate(subt);
      subt=0;
      la->block(b).assign_subblock(sub, nbf, 2*nbf-1, 0, nbf-1);

    } else {
      RefSymmSCMatrix irrep = la->block(b).copy();
      for (int c=1; c < ct.gamma(ir).degeneracy(); c++)
        irrep.accumulate(la->block(b+c));

      for (int c=0; c < ct.gamma(ir).degeneracy(); c++)
        la->block(b+c).assign(irrep);

      b += ct.gamma(ir).degeneracy()-1;
    }
  }
}


void
sc::symmetrize(const Ref<GPetiteList2>& plist12,
               const Ref<Integral>& integral,
               const RefSCMatrix& skel,
               const RefSCMatrix& sym)
{
  // validate input
  if (skel.rowdim().n() != plist12->basis1()->nbasis())
    throw ProgrammingError("sc::symmetrize(RefSCMatrix) -- row dimension of skel doesn't match the number of basis functions",__FILE__,__LINE__);
  if (skel.coldim().n() != plist12->basis2()->nbasis())
    throw ProgrammingError("sc::symmetrize(RefSCMatrix) -- col dimension of skel doesn't match the number of basis functions",__FILE__,__LINE__);
  if (skel.rowdim().n() != sym.rowdim().n())
    throw ProgrammingError("sc::symmetrize(RefSCMatrix) -- row dimensions of skel and sym don't match",__FILE__,__LINE__);
  if (skel.coldim().n() != sym.coldim().n())
    throw ProgrammingError("sc::symmetrize(RefSCMatrix) -- col dimensions of skel and sym don't match",__FILE__,__LINE__);

  const Ref<GaussianBasisSet>& bs1 = plist12->basis1();
  const Ref<GaussianBasisSet>& bs2 = plist12->basis2();
  integral->set_basis(bs1);
  Ref<PetiteList> pl1 = integral->petite_list();
  integral->set_basis(bs2);
  Ref<PetiteList> pl2 = integral->petite_list();

  // SO basis is always blocked, so first make sure skel is blocked
  RefSCMatrix bskel = dynamic_cast<BlockedSCMatrix*>(skel.pointer());
  if (bskel.null()) {
    bskel = bs1->so_matrixkit()->matrix(pl1->AO_basisdim(),pl2->AO_basisdim());
    bskel->convert(skel);
  }

  // if C1, then do nothing
  CharacterTable ct = plist12->point_group()->char_table();
  if (ct.order() == 1) {
    sym.assign(bskel);
    return;
  }

  RefSCMatrix aoso1 = pl1->aotoso();
  BlockedSCMatrix *lu1 = dynamic_cast<BlockedSCMatrix*>(aoso1.pointer());
  for (int b=0; b < lu1->nblocks(); b++) {
    if (lu1->block(b).null())
      continue;

    const int ir = ct.which_irrep(b);

    double skal = (double)ct.order()/(double)ct.gamma(ir).degeneracy();
    skal = sqrt(skal);
    lu1->block(b).scale(skal);
  }
  RefSCMatrix aoso2 = pl2->aotoso();
  BlockedSCMatrix *lu2 = dynamic_cast<BlockedSCMatrix*>(aoso2.pointer());
  for (int b=0; b < lu2->nblocks(); b++) {
    if (lu2->block(b).null())
      continue;

    const int ir = ct.which_irrep(b);

    double skal = (double)ct.order()/(double)ct.gamma(ir).degeneracy();
    skal = sqrt(skal);
    lu2->block(b).scale(skal);
  }

  sym.assign(0.0);
  RefSCMatrix tmp = aoso1.t() * bskel;
  sym.accumulate_product(tmp,aoso2);
  tmp = 0;
  aoso1 = 0;
  aoso2 = 0;

  // loop through blocks and finish symmetrizing degenerate blocks
  BlockedSCMatrix *la = dynamic_cast<BlockedSCMatrix*>(sym.pointer());
  for (int b=0; b < la->nblocks(); b++) {
    if (la->block(b).null())
      continue;

    const int ir=ct.which_irrep(b);

    if (ct.gamma(ir).degeneracy()==1)
      continue;

    if (ct.gamma(ir).complex()) {
      const int nbf1 = pl1->nfunction(ir)/2;
      const int nbf2 = pl2->nfunction(ir)/2;

      // average diagonal blocks
      RefSCMatrix irrep = la->block(b).get_subblock(0, nbf1-1, 0, nbf2-1);
      irrep.accumulate(la->block(b).get_subblock(nbf1, 2*nbf1-1, nbf2, 2*nbf2-1));
      la->block(b).assign_subblock(irrep,  0, nbf1-1, 0, nbf2-1);
      la->block(b).assign_subblock(irrep,nbf1, 2*nbf1-1, nbf2, 2*nbf2-1);

      // TODO figure out how to symmetrize rectangular matrices for the case of complex irreps
#if 0
      RefSCMatrix sub = la->block(b).get_subblock(nbf, 2*nbf-1, 0, nbf-1);
      RefSCMatrix subt = sub.t();
      subt.scale(-1.0);
      sub.accumulate(subt);
      subt=0;
      la->block(b).assign_subblock(sub, nbf, 2*nbf-1, 0, nbf-1);
#else
      abort();
#endif
    } else {
      RefSCMatrix irrep = la->block(b).copy();
      for (int c=1; c < ct.gamma(ir).degeneracy(); c++)
        irrep.accumulate(la->block(b+c));

      for (int c=0; c < ct.gamma(ir).degeneracy(); c++)
        la->block(b+c).assign(irrep);

      b += ct.gamma(ir).degeneracy()-1;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////

Ref<GPetiteList2>
sc::GPetiteListFactory::plist2(const Ref<GaussianBasisSet> &b1,
                               const Ref<GaussianBasisSet> &b2)
{
  if (b1 == b2) {
    return new GenericPetiteList2<canonical_aa>(b1,b2);
  }
  else {
    return new GenericPetiteList2<canonical_ab>(b1,b2);
  }
}

Ref<GPetiteList4>
sc::GPetiteListFactory::plist4(const Ref<GaussianBasisSet> &b1,
                  const Ref<GaussianBasisSet> &b2,
                  const Ref<GaussianBasisSet> &b3,
                  const Ref<GaussianBasisSet> &b4)
{
  if (b1 == b2 && b1 == b3 && b1 == b4) {
    return new GenericPetiteList4<canonical_aaaa>(b1,b2,b3,b4);
  }
  else if (b1 == b3 && b2 == b4) {
    return new GenericPetiteList4<canonical_abab>(b1,b2,b3,b4);
  }
  else if (b1 == b2 && b3 != b4) {
    return new GenericPetiteList4<canonical_aabc>(b1,b2,b3,b4);
  }
  else if (b1 != b2 && b3 == b4) {
    return new GenericPetiteList4<canonical_abcc>(b1,b2,b3,b4);
  }
  else if (b1 == b3 && b2 == b4) {
    return new GenericPetiteList4<canonical_abab>(b1,b2,b3,b4);
  }
  else if (b1 == b2 && b3 == b4) {
    return new GenericPetiteList4<canonical_aabb>(b1,b2,b3,b4);
  }
  else {
    return new GenericPetiteList4<canonical_abcd>(b1,b2,b3,b4);
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
