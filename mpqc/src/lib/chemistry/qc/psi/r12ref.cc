//
// r12ref.cc
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
#pragma implementation
#endif

#include <chemistry/qc/psi/r12ref.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/mbptr12/spin.h>

using namespace sc;

///////////////////////////////////////////////////////////////////

static ClassDesc PsiSD_R12RefWavefunction_cd(
  typeid(PsiSD_R12RefWavefunction),"PsiSD_R12RefWavefunction",1,"public R12RefWavefunction",
  0, 0, create<PsiSD_R12RefWavefunction>);

PsiSD_R12RefWavefunction::PsiSD_R12RefWavefunction(const Ref<WavefunctionWorld>& world,
                                                   const Ref<PsiSCF>& scf,
                                                   bool spin_restricted,
                                                   unsigned int nfzc,
                                                   unsigned int nfzv,
                                                   Ref<OrbitalSpace> vir_space) :
                                                   R12RefWavefunction(world,
                                                                      scf->basis(),
                                                                      scf->integral()),
                                                   scf_(scf),
                                                   spin_restricted_(spin_restricted),
                                                   nfzc_(nfzc),
                                                   nfzv_(nfzv),
                                                   vir_space_(vir_space) {
  // bring spin_restricted in sync with scf_
  if (scf_->spin_polarized() == false) spin_restricted_ = true;
  Ref<PsiUHF> uhf; uhf << scf_; if (uhf.nonnull()) spin_restricted_ = false;
  // make sure that FockBuildRuntime uses same densities as the reference wavefunction
  scf_->compute();
  world->fockbuild_runtime()->set_densities(this->ordm(Alpha), this->ordm(Beta));

  if (nfzv > 0 && vir_space.nonnull())
    throw ProgrammingError("when VBS is given nfzv must be 0",__FILE__,__LINE__);
}

PsiSD_R12RefWavefunction::PsiSD_R12RefWavefunction(StateIn& si) : R12RefWavefunction(si) {
  scf_ << SavableState::restore_state(si);
  vir_space_ << SavableState::restore_state(si);
  si.get(spin_restricted_);
  si.get(nfzc_);
  si.get(nfzv_);
}

PsiSD_R12RefWavefunction::~PsiSD_R12RefWavefunction() {
}

void
PsiSD_R12RefWavefunction::save_data_state(StateOut& so) {
  SavableState::save_state(scf_.pointer(), so);
  SavableState::save_state(vir_space_.pointer(), so);
  so.put(spin_restricted_);
  so.put(nfzc_);
  so.put(nfzv_);
}

namespace {
  SpinCase1 valid_spincase(SpinCase1 s) {
    return s == AnySpinCase1 ? Alpha : s;
  }
}

RefSymmSCMatrix
PsiSD_R12RefWavefunction::ordm(SpinCase1 s) const {
  s = valid_spincase(s);
  if (spin_restricted()) s = Alpha;
  RefSymmSCMatrix result = (s == Alpha) ? scf()->alpha_ao_density() : scf()->beta_ao_density();
  return result;
}

RefSymmSCMatrix
PsiSD_R12RefWavefunction::core_hamiltonian_for_basis(const Ref<GaussianBasisSet> &basis,
                                                  const Ref<GaussianBasisSet> &p_basis) {
  return scf()->core_hamiltonian_for_basis(basis, 0);
}

void
PsiSD_R12RefWavefunction::init_spaces()
{
  if (spin_restricted())
    init_spaces_restricted();
  else
    init_spaces_unrestricted();
}

void
PsiSD_R12RefWavefunction::init_spaces_restricted()
{
  const bool moorder = true;   // order orbitals in the order of increasing energies
  Ref<FockBuildRuntime> fbrun = this->world()->fockbuild_runtime();
  const Ref<GaussianBasisSet> bs = scf()->basis();
  const Ref<Integral>& integral =  scf()->integral();
  Ref<PetiteList> plist = integral->petite_list();
  RefSCMatrix evecs_ao = scf()->coefs();
  const RefDiagSCMatrix evals = scf()->evals();
  const int nmo = evecs_ao.coldim().n();

  // compute active orbital mask
  typedef MolecularOrbitalMask<double, RefDiagSCMatrix> FZCMask;
  typedef MolecularOrbitalMask<double, RefDiagSCMatrix, std::greater<double> > FZVMask;
  FZCMask fzcmask(nfzc(), evals);
  FZVMask fzvmask(nfzv(), evals);
  std::vector<bool> actmask(nmo, true);
  // add frozen core and frozen virtuals masks
  std::transform(fzcmask.mask().begin(), fzcmask.mask().end(),
                 fzvmask.mask().begin(), actmask.begin(), std::logical_and<bool>());

  using std::vector;
  vector<double> aoccs(nmo);
  for(int mo=0; mo<nmo; mo++) {
    aoccs[mo] = scf()->alpha_occupation(mo);
  }
  if (scf()->spin_polarized() == false) { // closed-shell
    spinspaces_[Alpha] = new PopulatedOrbitalSpace(AnySpinCase1, bs, integral, evecs_ao,
                                                   aoccs, actmask, evals, moorder,
                                                   vir_space(), fbrun);
    spinspaces_[Beta] = spinspaces_[Alpha];
  }
  else { // spin-restricted open-shell
    spinspaces_[Alpha] = new PopulatedOrbitalSpace(Alpha, bs, integral, evecs_ao,
                                                   aoccs, actmask, evals, moorder,
                                                   vir_space(), fbrun);
    vector<double> boccs(nmo);
    for(int mo=0; mo<nmo; mo++) {
      boccs[mo] = scf()->beta_occupation(mo);
    }
    spinspaces_[Beta] = new PopulatedOrbitalSpace(Beta, bs, integral, evecs_ao,
                                                  boccs, actmask, evals, moorder,
                                                  vir_space(), fbrun);
  }
}


void
PsiSD_R12RefWavefunction::init_spaces_unrestricted()
{
  const bool moorder = true;   // order orbitals in the order of increasing energies
  Ref<FockBuildRuntime> fbrun = this->world()->fockbuild_runtime();
  const Ref<GaussianBasisSet> bs = scf()->basis();
  const Ref<Integral>& integral = scf()->integral();

  using std::vector;
  vector<double> aocc, bocc;
  const int nmo = scf()->nmo();
  for(int mo=0; mo<nmo; mo++) {
    aocc.push_back(scf()->alpha_occupation(mo));
    bocc.push_back(scf()->beta_occupation(mo));
  }
  Ref<PetiteList> plist = scf()->integral()->petite_list();

  RefSCMatrix alpha_evecs, beta_evecs;
  RefDiagSCMatrix alpha_evals, beta_evals;
  // alpha and beta orbitals are available for UHF
  Ref<PsiUHF> uhf = dynamic_cast<PsiUHF*>(scf().pointer());
  if (uhf.nonnull()) {
    alpha_evecs = scf()->coefs(Alpha);
    beta_evecs = scf()->coefs(Beta);
    alpha_evals = scf()->evals(Alpha);
    beta_evals = scf()->evals(Beta);
  }
  // use semicanonical orbitals for ROHF
  else {
    Ref<PsiHSOSHF> hsoshf = dynamic_cast<PsiHSOSHF*>(ref().pointer());
    if (hsoshf.null())
      throw ProgrammingError("spin-specific spaces not available for this reference function", __FILE__, __LINE__);
    throw FeatureNotImplemented("semicanonical orbitals not yet implemented for PsiHSOSHF",__FILE__,__LINE__);
  }

  typedef MolecularOrbitalMask<double, RefDiagSCMatrix> FZCMask;
  typedef MolecularOrbitalMask<double, RefDiagSCMatrix, std::greater<double> > FZVMask;
  { // alpha spin
    // compute active orbital mask
    FZCMask fzcmask(nfzc(), alpha_evals);
    FZVMask fzvmask(nfzv(), alpha_evals);
    std::vector<bool> actmask(nmo, true);
    // add frozen core and frozen virtuals masks
    std::transform(fzcmask.mask().begin(), fzcmask.mask().end(),
                   fzvmask.mask().begin(), actmask.begin(), std::logical_and<bool>());
    spinspaces_[Alpha] = new PopulatedOrbitalSpace(Alpha, bs, integral,
                                                   alpha_evecs,
                                                   aocc, actmask, alpha_evals, moorder,
                                                   vir_space(), fbrun);
  }
  { // beta spin
    // compute active orbital mask
    FZCMask fzcmask(nfzc(), beta_evals);
    FZVMask fzvmask(nfzv(), beta_evals);
    std::vector<bool> actmask(nmo, true);
    // add frozen core and frozen virtuals masks
    std::transform(fzcmask.mask().begin(), fzcmask.mask().end(),
                   fzvmask.mask().begin(), actmask.begin(), std::logical_and<bool>());
    spinspaces_[Beta] = new PopulatedOrbitalSpace(Beta, bs, integral,
                                                  beta_evecs,
                                                  bocc, actmask, beta_evals, moorder,
                                                  vir_space(), fbrun);
  }
}

Ref<R12RefWavefunction>
make_PsiSD_R12RefWavefunction(const Ref<WavefunctionWorld>& world,
                              const Ref<Wavefunction>& wfn,
                              bool spin_restricted = true,
                              unsigned int nfzc = 0,
                              unsigned int nfzv = 0,
                              Ref<OrbitalSpace> vir_space = 0) {
  Ref<PsiSCF> scf; scf << wfn;
  assert(scf.nonnull());
  return new PsiSD_R12RefWavefunction(world, scf, spin_restricted, nfzc, nfzv, vir_space);
}


/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
