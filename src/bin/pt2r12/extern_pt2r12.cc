//
// extern_pt2r12.cc
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

#include <string>
#include "extern_pt2r12.h"
#include <iostream>

using namespace sc;

ClassDesc
ExternPT2R12::class_desc_(typeid(ExternPT2R12),
                     "ExternPT2R12",
                     1,               // version
                     "public Wavefunction",
                     0,               // this class is not DefaultConstructible
                     create<ExternPT2R12>, // this class is KeyValConstructible
                     0  // this class is not StateInConstructible
                     );

ExternPT2R12::ExternPT2R12(const Ref<KeyVal>& kv) :
    Wavefunction(kv)
{
  world_ << kv->describedclassvalue("world");
  world_->set_wfn(this);

  orbs_info_ << kv->describedclassvalue("orbs_info");
  rdm2_ << kv->describedclassvalue("rdm2");
  cabs_name_ = kv->stringvalue("cabs", KeyValValuestring(std::string()));
  f12exp_str_ = kv->stringvalue("f12exp", KeyValValuestring(std::string()));

  std::string r12_str = kv->stringvalue("pt2_correction", KeyValValuestring(std::string()));

#if defined(HAVE_MPQC3_RUNTIME)
  std::string mpqc3_str = kv->stringvalue("use_mpqc3", KeyValValuestring(std::string()));
  std::string singles_str = kv->stringvalue("cabs_singles", KeyValValuestring(std::string()));
  std::string partition_str = kv->stringvalue("cabs_singles_h0", KeyValValuestring(std::string()));
#endif

  Ref<OrbitalSpace> orbs = orbs_info_->orbs();
  const std::vector<unsigned int>& fzcpi = orbs_info_->fzcpi();
  const std::vector<unsigned int>& inactpi = orbs_info_->inactpi();
  const std::vector<unsigned int>& actpi = orbs_info_->actpi();
  const std::vector<unsigned int>& holepi = orbs_info_->corrpi();
  const std::vector<unsigned int>& fzvpi = orbs_info_->fzvpi();
  const std::vector<unsigned int>& mopi = orbs_info_->mopi();
  std::vector<unsigned int> occpi;
  const unsigned int nirrep = fzcpi.size();
  for (int i = 0; i < nirrep; ++i)
  {
    occpi.push_back(fzcpi[i] + inactpi[i] + actpi[i]);
  }

//  const unsigned int nfzc = std::accumulate(fzcpi.begin(), fzcpi.end(), 0.0);
//  const unsigned int ninact = std::accumulate(inactpi.begin(), inactpi.end(), 0.0);
//  const unsigned int nact = std::accumulate(actpi.begin(), actpi.end(), 0.0);
//  const unsigned int nfzv = std::accumulate(fzvpi.begin(), fzvpi.end(), 0.0);
//  const unsigned int nmo = orbs->rank();
//  const unsigned int nuocc = nmo - nfzc - ninact - nact;

  // to construct Extern_RefWavefunction we need an energy/correlation-ordered orbitals
  // make 1-RDM to the full MO space
  RefSymmSCMatrix P1_mo;
  {
    Ref<SpinFreeRDM<One> > rdrdm1 = rdm2_->rdm_m_1();
    RefSymmSCMatrix P1_mo_occ = rdrdm1->scmat();
    std::vector<unsigned int> occ_to_orbs_indexmap;
    occ_to_orbs_indexmap = (*orbs) << *(rdrdm1->orbs());
    P1_mo = P1_mo_occ.kit()->symmmatrix(orbs->dim());
    P1_mo.assign(0.0);
    const unsigned int nocc = P1_mo_occ.n();
    for(unsigned int i1=0; i1<nocc; ++i1) {
      const unsigned int ii1 = occ_to_orbs_indexmap[i1];
      for(unsigned int i2=0; i2<=i1; ++i2) {
        const unsigned int ii2 = occ_to_orbs_indexmap[i2];
        P1_mo.set_element(ii1, ii2, P1_mo_occ.get_element(i1, i2));
      }
    }
    P1_mo.scale(0.5);
  }

  // use its orbitals to initialize Extern_RefWavefunction
  Ref<Integral> intf = this->integral()->clone();
  intf->set_basis(basis());
  Ref<RefWavefunction> ref_wfn = new Extern_RefWavefunction(world_, basis(), intf,
                                                            orbs->coefs(), orbs->orbsym(),
                                                            P1_mo, P1_mo,
                                                            occpi,
                                                            fzcpi,
                                                            fzvpi,
                                                            holepi);
  if(debug_print_)
  {
    sc::ExEnv::out0() << "debug:print refwfn orbs " << std::endl;
    ref_wfn->occ_act_sb(Alpha)->print_detail();
    ref_wfn->orbs_sb(Alpha)->print_detail();
  }

  // create PT2R12 object
  {
    Ref<AssignedKeyVal> kva = new AssignedKeyVal;
    kva->assign("molecule", molecule().pointer());
    kva->assign("basis", basis().pointer());
    kva->assign("refwfn", ref_wfn.pointer());
    kva->assign("world", world_.pointer());
    kva->assign("rdm2", rdm2_.pointer());
    kva->assign("corr_factor", "stg-6g");
    kva->assign("corr_param", f12exp_str_.c_str());
    if(!r12_str.empty())
      kva->assign("pt2_correction", r12_str);

#if defined(HAVE_MPQC3_RUNTIME)
    if(!singles_str.empty())
      kva->assign("cabs_singles", singles_str);
    if(!partition_str.empty())
      kva->assign("cabs_singles_h0", partition_str);
    if(!mpqc3_str.empty())
      kva->assign("use_mpqc3", mpqc3_str);
#endif
    if (cabs_name_.empty() == false) {
      Ref<AssignedKeyVal> tmpkv = new AssignedKeyVal;
      tmpkv->assign("name", cabs_name_.c_str());
      tmpkv->assign("puream", "true");
      tmpkv->assign("molecule", molecule().pointer());
      Ref<KeyVal> kv = tmpkv;
      Ref<GaussianBasisSet> aux_basis = new GaussianBasisSet(kv);
      kva->assign("aux_basis", aux_basis.pointer());
    }
    else { // CABS name not given? construct automatically using R12Technology::make_auto_cabs()
      kva->assign("aux_basis", R12Technology::make_auto_cabs(basis()).pointer());
    }
    Ref<KeyVal> kv = kva;

    pt2r12_ = new PT2R12(kv);
  }
}

ExternPT2R12::~ExternPT2R12() {
  const bool make_sure_class_desc_initialized = (&class_desc_ != 0);
}

void ExternPT2R12::compute()
{
  const double value = pt2r12_->value();
}

int ExternPT2R12::nelectron()
{
  return pt2r12_->nelectron();
}

RefSymmSCMatrix ExternPT2R12::density()
{
  throw FeatureNotImplemented("density evaluation for ExternPT2R12 not yet implemented",
                              __FILE__, __LINE__, this->class_desc());
}

double ExternPT2R12::magnetic_moment() const
{
  return pt2r12_->magnetic_moment();
}

void ExternPT2R12::set_desired_value_accuracy(double acc) {
  pt2r12_->set_desired_value_accuracy(acc);
}

void ExternPT2R12::print(std::ostream& os) {
  pt2r12_->print(os);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
