//
// singlerefinfo.cc
//
// Copyright (C) 2005 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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

#include <util/misc/scexception.h>
#include <chemistry/qc/scf/uhf.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/mbptr12/singlerefinfo.h>

using namespace sc;

namespace {
  const unsigned int ClassVersion = 1;
};

static ClassDesc R12IntEvalInfo_cd(
  typeid(SingleRefInfo),"SingleRefInfo",ClassVersion,"virtual public SavableState",
  0, 0, create<SingleRefInfo>);

SingleRefInfo::SingleRefInfo(const Ref<SCF>& ref, unsigned int nfzc, unsigned int nfzv) :
  ref_(ref), nfzc_(nfzc), nfzv_(nfzv)
{
  if (!ref_->spin_polarized())
    init_spinindependent_spaces();
  init_spinspecific_spaces();
}

SingleRefInfo::SingleRefInfo(StateIn& si) :
  SavableState(si)
{
  ref_ << SavableState::restore_state(si);
  si.get(nfzc_);
  si.get(nfzv_);
  
  orbs_sb_ << SavableState::restore_state(si);
  orbs_ << SavableState::restore_state(si);
  docc_sb_ << SavableState::restore_state(si);
  docc_ << SavableState::restore_state(si);
  docc_act_ << SavableState::restore_state(si);
  socc_ << SavableState::restore_state(si);
  uocc_sb_ << SavableState::restore_state(si);
  uocc_ << SavableState::restore_state(si);
  uocc_act_ << SavableState::restore_state(si);
  for(int spin=0; spin<2; spin++) {
    spinspaces_[spin].orbs_sb_ << SavableState::restore_state(si);
    spinspaces_[spin].orbs_ << SavableState::restore_state(si);
    spinspaces_[spin].occ_sb_ << SavableState::restore_state(si);
    spinspaces_[spin].occ_ << SavableState::restore_state(si);
    spinspaces_[spin].occ_act_ << SavableState::restore_state(si);
    spinspaces_[spin].uocc_sb_ << SavableState::restore_state(si);
    spinspaces_[spin].uocc_ << SavableState::restore_state(si);
    spinspaces_[spin].uocc_act_ << SavableState::restore_state(si);
  }
}

SingleRefInfo::~SingleRefInfo()
{
}

void
SingleRefInfo::save_data_state(StateOut& so)
{
  SavableState::save_state(ref_.pointer(),so);
  so.put(nfzc_);
  so.put(nfzv_);
  
  SavableState::save_state(orbs_sb_.pointer(),so);
  SavableState::save_state(orbs_.pointer(),so);
  SavableState::save_state(docc_sb_.pointer(),so);
  SavableState::save_state(docc_.pointer(),so);
  SavableState::save_state(docc_act_.pointer(),so);
  SavableState::save_state(socc_.pointer(),so);
  SavableState::save_state(uocc_sb_.pointer(),so);
  SavableState::save_state(uocc_.pointer(),so);
  SavableState::save_state(uocc_act_.pointer(),so);
  for(int spin=0; spin<2; spin++) {
    SavableState::save_state(spinspaces_[spin].orbs_sb_.pointer(),so);
    SavableState::save_state(spinspaces_[spin].orbs_.pointer(),so);
    SavableState::save_state(spinspaces_[spin].occ_sb_.pointer(),so);
    SavableState::save_state(spinspaces_[spin].occ_.pointer(),so);
    SavableState::save_state(spinspaces_[spin].occ_act_.pointer(),so);
    SavableState::save_state(spinspaces_[spin].uocc_sb_.pointer(),so);
    SavableState::save_state(spinspaces_[spin].uocc_.pointer(),so);
    SavableState::save_state(spinspaces_[spin].uocc_act_.pointer(),so);
  }
}

const Ref<SCF>&
SingleRefInfo::ref() const
{
  return ref_;
}

unsigned int
SingleRefInfo::nfzc() const
{
  return nfzc_;
}

unsigned int
SingleRefInfo::nfzv() const
{
  return nfzv_;
}

void
SingleRefInfo::init_spinspecific_spaces()
{
  const Ref<GaussianBasisSet> bs = ref_->basis();
  using std::vector;
  vector<double> aocc, bocc;
  const int nmo = ref_->alpha_eigenvectors().coldim().n();
  for(int mo=0; mo<nmo; mo++) {
    aocc.push_back(ref_->alpha_occupation(mo));
    bocc.push_back(ref_->beta_occupation(mo));
  }
  Ref<PetiteList> plist = ref_->integral()->petite_list();
  if (ref()->spin_polarized()) {
    spinspaces_[0].init("Alpha", bs, ref_->alpha_eigenvalues(), plist->evecs_to_AO_basis(ref_->alpha_eigenvectors()), aocc, nfzc(), nfzv());
    spinspaces_[1].init("Beta", bs, ref_->beta_eigenvalues(), plist->evecs_to_AO_basis(ref_->beta_eigenvectors()), bocc, nfzc(), nfzv());
  }
  else {
    for(int s=0; s<NSpinCases1; s++) {
      spinspaces_[s].orbs_sb_ = orbs_sb_;
      spinspaces_[s].orbs_ = orbs_;
      spinspaces_[s].occ_sb_ = docc_sb();
      spinspaces_[s].occ_ = docc();
      spinspaces_[s].occ_act_ = docc_act();
      spinspaces_[s].uocc_sb_ = uocc_sb_;
      spinspaces_[s].uocc_ = uocc_;
      spinspaces_[s].uocc_act_ = uocc_act_;
    }
  }
}

void
SingleRefInfo::init_spinindependent_spaces()
{
  const Ref<GaussianBasisSet> bs = ref_->basis();
  const RefSCMatrix evecs_so = ref_->eigenvectors();
  const RefDiagSCMatrix evals = ref_->eigenvalues();
  Ref<PetiteList> plist = ref_->integral()->petite_list();
  RefSCMatrix evecs_ao = plist->evecs_to_AO_basis(evecs_so);
  
  int ndocc = 0;
  int nsocc = 0;
  int nuocc = 0;
  const int nmo = evecs_ao.coldim().n();
  for (int i=0; i<nmo; i++) {
    if (ref_->occupation(i) == 2.0)
      ndocc++;
    else if (ref_->occupation(i) == 1.0)
      nsocc++;
    else
      nuocc++;
  }
  
  orbs_sb_ = new MOIndexSpace("p(sym)","symmetry-blocked MOs", evecs_ao, bs, evals, 0, 0, MOIndexSpace::symmetry);
  orbs_ = new MOIndexSpace("p","energy-ordered MOs", evecs_ao, bs, evals, 0, 0);
  docc_sb_ = new MOIndexSpace("m(sym)","doubly-occupied symmetry-blocked MOs", evecs_ao, bs, evals, 0, nuocc+nsocc, MOIndexSpace::symmetry);
  docc_ = new MOIndexSpace("m","doubly-occupied energy-ordered MOs", evecs_ao, bs, evals, 0, nuocc+nsocc);
  docc_act_ = new MOIndexSpace("i","active doubly-occupied energy-ordered MOs", evecs_ao, bs, evals, nfzc(), nuocc+nsocc);
  socc_ = new MOIndexSpace("x","singly-occupied energy-ordered MOs", evecs_ao, bs, evals, ndocc, nuocc);
  uocc_sb_ = new MOIndexSpace("e(sym)","unoccupied symmetry-blocked MOs", evecs_ao, bs, evals, ndocc+nsocc, 0, MOIndexSpace::symmetry);
  uocc_ = new MOIndexSpace("e","unoccupied energy-ordered MOs", evecs_ao, bs, evals, ndocc+nsocc, 0);
  uocc_act_ = new MOIndexSpace("a","active unoccupied energy-ordered MOs", evecs_ao, bs, evals, ndocc+nsocc, nfzv());
}


const Ref<MOIndexSpace>&
SingleRefInfo::orbs_sb() const
{
  // can throw
  throw_if_spin_polarized();
  return orbs_sb_;
}

const Ref<MOIndexSpace>&
SingleRefInfo::orbs() const
{
  // can throw
  throw_if_spin_polarized();
  return orbs_;
}

const Ref<MOIndexSpace>&
SingleRefInfo::docc_sb() const
{
  // can throw
  throw_if_spin_polarized();
  return docc_sb_;
}

const Ref<MOIndexSpace>&
SingleRefInfo::docc() const
{
  // can throw
  throw_if_spin_polarized();
  return docc_;
}

const Ref<MOIndexSpace>&
SingleRefInfo::docc_act() const
{
  // can throw
  throw_if_spin_polarized();
  return docc_act_;
}
const Ref<MOIndexSpace>&
SingleRefInfo::socc() const
{
  // can throw
  throw_if_spin_polarized();
  return socc_;
}

const Ref<MOIndexSpace>&
SingleRefInfo::uocc_sb() const
{
  // can throw
  throw_if_spin_polarized();
  return uocc_sb_;
}

const Ref<MOIndexSpace>&
SingleRefInfo::uocc() const
{
  // can throw
  throw_if_spin_polarized();
  return uocc_;
}

const Ref<MOIndexSpace>&
SingleRefInfo::uocc_act() const
{
  // can throw
  throw_if_spin_polarized();
  return uocc_act_;
}

const Ref<MOIndexSpace>&
SingleRefInfo::orbs_sb(SpinCase1 s) const
{
  return spinspaces_[s].orbs_sb_;
}

const Ref<MOIndexSpace>&
SingleRefInfo::orbs(SpinCase1 s) const
{
  return spinspaces_[s].orbs_;
}

const Ref<MOIndexSpace>&
SingleRefInfo::occ_sb(SpinCase1 s) const
{
  return spinspaces_[s].occ_sb_;
}

const Ref<MOIndexSpace>&
SingleRefInfo::occ(SpinCase1 s) const
{
  return spinspaces_[s].occ_;
}

const Ref<MOIndexSpace>&
SingleRefInfo::occ_act(SpinCase1 s) const
{
  return spinspaces_[s].occ_act_;
}

const Ref<MOIndexSpace>&
SingleRefInfo::uocc_sb(SpinCase1 s) const
{
  return spinspaces_[s].uocc_sb_;
}

const Ref<MOIndexSpace>&
SingleRefInfo::uocc(SpinCase1 s) const
{
  return spinspaces_[s].uocc_;
}

const Ref<MOIndexSpace>&
SingleRefInfo::uocc_act(SpinCase1 s) const
{
  return spinspaces_[s].uocc_act_;
}

void
SingleRefInfo::throw_if_spin_polarized() const
{
  if (ref_->spin_polarized())
    throw ProgrammingError("SingleRefInfo -- spin-independent space is requested but the reference function is spin-polarized",
        __FILE__,__LINE__);
}

/////////////

void
SingleRefInfo::SpinSpaces::init(
        const std::string& prefix,
        const Ref<GaussianBasisSet>& bs,
        const RefDiagSCMatrix& evals,
        const RefSCMatrix& evecs,
        const std::vector<double>& occs,
        unsigned int nfzc,
        unsigned int nfzv)
{
  int nocc = 0, nuocc = 0, nmo = occs.size();
  for(int i=0; i<nmo; i++) {
    if (occs[i] == 1.0)
      nocc++;
    else
      nuocc++;
  }
  using std::ostringstream;
  {
    ostringstream oss;
    oss << prefix << " symmetry-blocked MOs";
    std::string id = (prefix[0]=='A'?"P(sym)":"p(sym)");
    orbs_sb_ = new MOIndexSpace(id,oss.str(),evecs, bs, evals, 0, 0, MOIndexSpace::symmetry);
  }
  {
    ostringstream oss;
    oss << prefix << " energy-ordered MOs";
    std::string id = (prefix[0]=='A'?"P":"p");
    orbs_ = new MOIndexSpace(id,oss.str(),evecs, bs, evals, 0, 0);
  }
  {
    ostringstream oss;
    oss << prefix << " occupied symmetry-blocked MOs";
    std::string id = (prefix[0]=='A'?"M(sym)":"m(sym)");
    occ_sb_ = new MOIndexSpace(id,oss.str(),evecs, bs, evals, 0, nuocc, MOIndexSpace::symmetry);
  }
  {
    ostringstream oss;
    oss << prefix << " occupied MOs";
    std::string id = (prefix[0]=='A'?"M":"m");
    occ_ = new MOIndexSpace(id,oss.str(),evecs, bs, evals, 0, nuocc);
  }
  {
    ostringstream oss;
    oss << prefix << " active occupied MOs";
    std::string id = (prefix[0]=='A'?"I":"i");
    occ_act_ = new MOIndexSpace(id,oss.str(),evecs, bs, evals, nfzc, nuocc);
  }
  {
    ostringstream oss;
    oss << prefix << " unoccupied symmetry-blocked MOs";
    std::string id = (prefix[0]=='A'?"E(sym)":"e(sym)");
    uocc_sb_ = new MOIndexSpace(id,oss.str(),evecs, bs, evals, nocc, 0, MOIndexSpace::symmetry);
  }
  {
    ostringstream oss;
    oss << prefix << " unoccupied MOs";
    std::string id = (prefix[0]=='A'?"E":"e");
    uocc_ = new MOIndexSpace(id,oss.str(),evecs, bs, evals, nocc, 0);
  }
  {
    ostringstream oss;
    oss << prefix << " active unoccupied MOs";
    std::string id = (prefix[0]=='A'?"A":"a");
    uocc_act_ = new MOIndexSpace(id,oss.str(),evecs, bs, evals, nocc, nfzv);
  }
}


