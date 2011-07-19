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
    tmpkv->assign("name", "mod6-31Gs");
    tmpkv->assign("basisdir", "./");
    tmpkv->assign("molecule", molecule.pointer());
    Ref<KeyVal> kv = tmpkv;
    Ref<GaussianBasisSet> basis = new GaussianBasisSet(kv);

#if 0
    // split generally contracted shells -- IntegralLibint2 can't handle such shells
    Ref<AssignedKeyVal> tmpkv1 = new AssignedKeyVal;
    tmpkv1->assign("basis", basis.pointer());
    kv = tmpkv1;
    basis_ = new SplitBasisSet(kv);
#endif
    basis_ = basis;
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

  RefSCDimension aodim = new SCDimension(nao, 1);
  aodim->blocks()->set_subdim(0, new SCDimension(nao));
  RefSCDimension modim = new SCDimension(nmo, 1);
  modim->blocks()->set_subdim(0, new SCDimension(nmo));
  coefs_ = basis_->so_matrixkit()->matrix(aodim, modim);
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
  if (in.is_open() == false) {
    std::ostringstream oss;
    oss << "ExternReadRDMOne: could not open file " << filename;
    throw std::runtime_error(oss.str().c_str());
  }
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
  in.close();

  rdm_.print("ExternReadRDMOne:: MO density");
}

ExternReadRDMOne::ExternReadRDMOne(const RefSymmSCMatrix& rdm, const Ref<OrbitalSpace>& orbs) :
  RDM<One>(Ref<Wavefunction>()), rdm_(rdm), orbs_(orbs)
{
  rdm_.print("ExternReadRDMOne:: MO density");
}

ExternReadRDMOne::~ExternReadRDMOne()
{
}

/////////////////////////////////////////////////////////////////////////////

ClassDesc ExternReadRDMTwo::class_desc_(typeid(ExternReadRDMTwo),
                                        "ExternReadRDMTwo",
                                        1,               // version
                                        "public RDM<Two>", // must match parent
                                        0,               // not DefaultConstructible
                                        0,               // not KeyValConstructible
                                        0                // not StateInConstructible
                                        );

sc::ExternReadRDMTwo::ExternReadRDMTwo(const std::string & filename, const Ref<OrbitalSpace> & orbs) :
    RDM<Two>(Ref<Wavefunction>()), filename_(filename), orbs_(orbs)
{
  std::ifstream in(filename_.c_str());
  if (in.is_open() == false) {
    std::ostringstream oss;
    oss << "ExternReadRDMTwo: could not open file " << filename_;
    throw std::runtime_error(oss.str().c_str());
  }
  const unsigned int norbs = orbs->coefs().coldim().n();
  RefSCDimension dim = new SCDimension(norbs * norbs);
  dim->blocks()->set_subdim(0, new SCDimension(dim->n()));
  rdm_ = orbs->coefs().kit()->symmmatrix(dim);
  rdm_.assign(0.0);
  bool have_coefs = true;
  while (have_coefs) {
    int bra1, bra2, ket1, ket2;
    double value;
    in >> bra1 >> ket1 >> ket2 >> bra2 >> value;
    if (bra1 != -1) {
      --bra1;
      --bra2;
      --ket1;
      --ket2;
      
      // this is specific to Luke's early RDM2 files produced from GAMESS
      value *= 4.0;
      if (bra1 == bra2 && ket1 == ket2 && bra1 == ket1) value *= 2.0;
      //
      
      rdm_.set_element(bra1 * norbs + bra2, ket1 * norbs + ket2, value);
      rdm_.set_element(bra2 * norbs + bra1, ket2 * norbs + ket1, value);
    } else
      have_coefs = false;
  }
  in.close();

  //rdm_.print("ExternReadRDMTwo:: MO density");
}

ExternReadRDMTwo::~ExternReadRDMTwo()
{
}

Ref<ExternReadRDMTwo::cumulant_type>
ExternReadRDMTwo::cumulant() const
{
  return new RDMCumulant<Two>(const_cast<ExternReadRDMTwo*>(this));
}

Ref< RDM<One> >
ExternReadRDMTwo::rdm_m_1() const
{
  Ref<ExternReadRDMOne> result;
  bool have_rdm1_file = true;
  try {
    std::string rdm1_filename(filename_);
    // compute filename for rdm1
    // hardwire for now
    rdm1_filename = "crap.txt";
    result = new ExternReadRDMOne(rdm1_filename, orbs_);
  }
  catch (...) {
    have_rdm1_file = false;
  }
  // if could not find the rdm1 in a file, compute it
  if (have_rdm1_file == false) {
    RefSymmSCMatrix rdm1 = orbs_->coefs().kit()->symmmatrix(orbs_->coefs().coldim());
    rdm1.assign(0.0);
    const unsigned int norbs = rdm1.n();
    for(unsigned int b1=0; b1<norbs; ++b1) {
      const unsigned b12_offset = b1 * norbs;
      for(unsigned int k1=0; k1<=b1; ++k1) {
        const unsigned k12_offset = k1 * norbs;
        double value = 0.0;
        for(unsigned int i2=0; i2<norbs; ++i2) {
          value += rdm_.get_element(b12_offset + i2, k12_offset + i2);
        }
        rdm1.set_element(b1, k1, value);
      }
    }

    // compute the number of electrons and scale 1-rdm:
    // trace of the 2-rdm is n(n-1)
    // trace of the 1-rdm should be n
    const double trace_2rdm = rdm1.trace();
    const double nelectron = (1.0 + std::sqrt(1.0 + 4.0 * trace_2rdm)) / 2.0;
    rdm1.scale(1.0 / (nelectron - 1.0));

    result = new ExternReadRDMOne(rdm1, orbs_);
  }

  return result;
}





/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
