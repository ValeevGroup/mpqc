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

SingleRefInfo::SingleRefInfo(const Ref<SCF>& ref) :
  ref_(ref)
{
  if (!ref_->spin_polarized())
    init_spinindependent_spaces();
  init_spinspecific_spaces();
}

SingleRefInfo::SingleRefInfo(StateIn& si) :
  SavableState(si)
{
  ref_ << SavableState::restore_state(si);
  
  symblk_mo_ << SavableState::restore_state(si);
  energy_mo_ << SavableState::restore_state(si);
  docc_ << SavableState::restore_state(si);
  socc_ << SavableState::restore_state(si);
  uocc_ << SavableState::restore_state(si);
  for(int spin=0; spin<2; spin++) {
    spinspaces_[spin].symblk_mo_ << SavableState::restore_state(si);
    spinspaces_[spin].energy_mo_ << SavableState::restore_state(si);
    spinspaces_[spin].occ_ << SavableState::restore_state(si);
    spinspaces_[spin].uocc_ << SavableState::restore_state(si);
  }
}

void
SingleRefInfo::save_data_state(StateOut& so)
{
  SavableState::save_state(ref_.pointer(),so);
  
  SavableState::save_state(symblk_mo_.pointer(),so);
  SavableState::save_state(energy_mo_.pointer(),so);
  SavableState::save_state(docc_.pointer(),so);
  SavableState::save_state(socc_.pointer(),so);
  SavableState::save_state(uocc_.pointer(),so);
  for(int spin=0; spin<2; spin++) {
    SavableState::save_state(spinspaces_ [spin].symblk_mo_.pointer(),so);
    SavableState::save_state(spinspaces_[spin].energy_mo_.pointer(),so);
    SavableState::save_state(spinspaces_[spin].occ_.pointer(),so);
    SavableState::save_state(spinspaces_[spin].uocc_.pointer(),so);
  }
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
  spinspaces_[0].init("Alpha", bs,ref_->alpha_eigenvalues(),ref_->alpha_eigenvectors(),aocc);
  spinspaces_[1].init("Beta", bs,ref_->beta_eigenvalues(),ref_->beta_eigenvectors(),bocc);
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
  
  symblk_mo_ = new MOIndexSpace("symmetry-blocked MOs", evecs_ao, bs, evals, 0, 0, MOIndexSpace::symmetry);
  energy_mo_ = new MOIndexSpace("energy-ordered MOs", evecs_ao, bs, evals, 0, 0);
  docc_ = new MOIndexSpace("doubly-occupied energy-ordered MOs", evecs_ao, bs, evals, 0, nuocc+nsocc);
  socc_ = new MOIndexSpace("singly-occupied energy-ordered MOs", evecs_ao, bs, evals, ndocc, nuocc);
  uocc_ = new MOIndexSpace("unoccupied energy-ordered MOs", evecs_ao, bs, evals, ndocc+nsocc, 0);
  
}


const Ref<MOIndexSpace>&
SingleRefInfo::symblk_mo() const
{
  // can throw
  throw_if_spin_polarized();
  return symblk_mo_;
}

const Ref<MOIndexSpace>&
SingleRefInfo::energy_mo() const
{
  // can throw
  throw_if_spin_polarized();
  return energy_mo_;
}

const Ref<MOIndexSpace>&
SingleRefInfo::docc() const
{
  // can throw
  throw_if_spin_polarized();
  return docc_;
}

const Ref<MOIndexSpace>&
SingleRefInfo::socc() const
{
  // can throw
  throw_if_spin_polarized();
  return socc_;
}

const Ref<MOIndexSpace>&
SingleRefInfo::uocc() const
{
  // can throw
  throw_if_spin_polarized();
  return uocc_;
}

const Ref<MOIndexSpace>&
SingleRefInfo::symblk_mo(SpinCase s) const
{
  return spinspaces_[s].symblk_mo_;
}

const Ref<MOIndexSpace>&
SingleRefInfo::energy_mo(SpinCase s) const
{
  return spinspaces_[s].energy_mo_;
}

const Ref<MOIndexSpace>&
SingleRefInfo::occ(SpinCase s) const
{
  return spinspaces_[s].occ_;
}

const Ref<MOIndexSpace>&
SingleRefInfo::uocc(SpinCase s) const
{
  return spinspaces_[s].uocc_;
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
        const std::vector<double>& occs)
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
    symblk_mo_ = new MOIndexSpace(oss.str(),evecs, bs, evals, 0, 0, MOIndexSpace::symmetry);
  }
  {
    ostringstream oss;
    oss << prefix << " energy-ordered MOs";
    energy_mo_ = new MOIndexSpace(oss.str(),evecs, bs, evals, 0, 0);
  }
  {
    ostringstream oss;
    oss << prefix << " occupied MOs";
    occ_ = new MOIndexSpace(oss.str(),evecs, bs, evals, 0, nuocc);
  }
  {
    ostringstream oss;
    oss << prefix << " unoccupied MOs";
    uocc_ = new MOIndexSpace(oss.str(),evecs, bs, evals, nocc, 0);
  }
}


