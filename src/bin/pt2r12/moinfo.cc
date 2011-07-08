//
// moinfo.cc
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

#ifdef __GNUG__
#pragma implementation
#endif

#include <moinfo.h>
#include <iostream>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/basis/split.h>

using namespace sc;

ExternReadMOInfo::ExternReadMOInfo(const std::string & filename)
{
  std::ifstream in(filename.c_str());

  Ref<Molecule> molecule = new Molecule;
  bool have_atoms = true;
  while (have_atoms) {
    double charge, x, y, z;
    in >> charge >> x >> y >> z;
    if (charge != -1.0) {
      molecule->add_atom(charge, x, y, z);
    }
    else
      have_atoms = false;
  }
  //molecule->set_point_group(new PointGroup("d2h"));
  molecule->print();

  // initialize basis set given a hardwired name
  {
    Ref<AssignedKeyVal> tmpkv = new AssignedKeyVal;
    tmpkv->assign("name", "6-31G*");
    //tmpkv->assign("basisdir", "./");
    tmpkv->assign("molecule", molecule.pointer());
    Ref<KeyVal> kv = tmpkv;
    Ref<GaussianBasisSet> basis = new GaussianBasisSet(kv);

    // split generally contracted shells -- IntegralLibint2 can't handle such shells
    Ref<AssignedKeyVal> tmpkv1 = new AssignedKeyVal;
    tmpkv1->assign("basis", basis.pointer());
    kv = tmpkv1;
    basis_ = new SplitBasisSet(kv);
  }
  basis_->print();

  unsigned int nmo;  in >> nmo;
  unsigned int nao;  in >> nao;
  unsigned int ncore; in >> ncore;
  in >> nfzc_;
  unsigned int nact; in >> nact;
  nocc_ = nact + ncore;
  // for now assume no frozen virtuals
  nfzv_ = 0;

  orbsym_.resize(nmo);
  for(int o=0; o<nmo; ++o) {
    unsigned int irrep; in >> irrep;
    --irrep;
    orbsym_[o] = irrep;
  }
  unsigned int junk; in >> junk;

  Ref<Integral> integral = Integral::get_default_integral()->clone();
  integral->set_basis(basis_);
  Ref<PetiteList> pl = integral->petite_list();
  coefs_ = basis_->so_matrixkit()->matrix(pl->AO_basisdim(), pl->SO_basisdim());
  coefs_.assign(0.0);
  bool have_coefs = true;
  while(have_coefs) {
    int row, col;
    double value;
    in >> row >> col >> value;
    if (row != -1) {
      --row;  --col;
      coefs_.set_element(col,row,value);
    }
    else
      have_coefs = false;
  }

  coefs_.print("ExternReadMOInfo:: MO coefficients");

  in.close();
}


Ref<GaussianBasisSet> ExternReadMOInfo::basis() const
{
    return basis_;
}

RefSCMatrix ExternReadMOInfo::coefs() const
{
    return coefs_;
}

unsigned int ExternReadMOInfo::nfzc() const
{
    return nfzc_;
}

unsigned int ExternReadMOInfo::nfzv() const
{
    return nfzv_;
}

unsigned int ExternReadMOInfo::nocc() const
{
    return nocc_;
}

std::vector<unsigned int> ExternReadMOInfo::orbsym() const
{
    return orbsym_;
}

/////////////////////////////////////////////////////////////////////////////

ClassDesc ExternReadRDMOne::class_desc_(typeid(ExternReadRDMOne),
                                        "ExternReadRDMOne",
                                        1,               // version
                                        "public RDM<One>", // must match parent
                                        0,               // not DefaultConstructible
                                        0,               // not KeyValConstructible
                                        0                // not StateInConstructible
                                        );

ExternReadRDMOne::ExternReadRDMOne(const std::string & filename, const Ref<OrbitalSpace>& orbs) :
  RDM<One>(Ref<Wavefunction>()), orbs_(orbs)
{
  std::ifstream in(filename.c_str());

  rdm_ = orbs->coefs().kit()->symmmatrix(orbs->coefs().coldim()); rdm_.assign(0.0);
  bool have_coefs = true;
  while(have_coefs) {
    int row, col;
    double value;
    in >> row >> col >> value;
    if (row != -1) {
      --row;  --col;
      rdm_.set_element(row,col,value);
    }
    else
      have_coefs = false;
  }

  rdm_.print("ExternReadRDMOne:: MO density");
}

ExternReadRDMOne::~ExternReadRDMOne()
{
}



/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
