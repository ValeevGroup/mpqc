//
// psiref.cc
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

#include <chemistry/qc/psi/psiref.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/mbptr12/spin.h>

using namespace sc;

///////////////////////////////////////////////////////////////////

static ClassDesc PsiSCF_R12RefWavefunction_cd(
  typeid(PsiSCF_RefWavefunction),"PsiSCF_R12RefWavefunction",1,"public R12RefWavefunction",
  0, 0, create<PsiSCF_RefWavefunction>);

PsiSCF_RefWavefunction::PsiSCF_RefWavefunction(const Ref<WavefunctionWorld>& world,
                                                   const Ref<PsiSCF>& scf,
                                                   bool spin_restricted,
                                                   unsigned int nfzc,
                                                   unsigned int nfzv,
                                                   Ref<OrbitalSpace> vir_space) :
                                                   RefWavefunction(world,
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

  if (nfzv > 0 && vir_space.nonnull())
    throw ProgrammingError("when VBS is given nfzv must be 0",__FILE__,__LINE__);
}

PsiSCF_RefWavefunction::PsiSCF_RefWavefunction(StateIn& si) : RefWavefunction(si) {
  scf_ << SavableState::restore_state(si);
  vir_space_ << SavableState::restore_state(si);
  si.get(spin_restricted_);
  si.get(nfzc_);
  si.get(nfzv_);
}

PsiSCF_RefWavefunction::~PsiSCF_RefWavefunction() {
}

void
PsiSCF_RefWavefunction::save_data_state(StateOut& so) {
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
PsiSCF_RefWavefunction::ordm(SpinCase1 s) const {
  s = valid_spincase(s);
  if (!scf_->spin_polarized()) s = Alpha;
  RefSymmSCMatrix result = (s == Alpha) ? scf()->alpha_ao_density() : scf()->beta_ao_density();
  return result;
}

RefSymmSCMatrix
PsiSCF_RefWavefunction::core_hamiltonian_for_basis(const Ref<GaussianBasisSet> &basis,
                                                  const Ref<GaussianBasisSet> &p_basis) {
  return scf()->core_hamiltonian_for_basis(basis, 0);
}

void
PsiSCF_RefWavefunction::init_spaces()
{
  if (spin_restricted())
    init_spaces_restricted();
  else
    init_spaces_unrestricted();
}

void
PsiSCF_RefWavefunction::init_spaces_restricted()
{
  const bool moorder = true;   // order orbitals in the order of increasing energies
  Ref<FockBuildRuntime> fbrun = this->world()->fockbuild_runtime();
  const Ref<GaussianBasisSet> bs = scf()->basis();
  const Ref<Integral>& integral =  scf()->integral();
  Ref<PetiteList> plist = integral->petite_list();
  RefSCMatrix evecs_ao = scf()->coefs();
  RefDiagSCMatrix evals = scf()->evals();
  int nmo = evals.n();

  using std::vector;
  vector<double> aoccs(nmo);
  vector<double> boccs(nmo);
  for(int mo=0; mo<nmo; mo++) {
    aoccs[mo] = scf()->alpha_occupation(mo);
    boccs[mo] = scf()->beta_occupation(mo);
  }

  // omit unoccupied orbitals?
  const bool omit_uocc = vir_space_.nonnull() && vir_space_->rank() == 0;
  if (omit_uocc) {
    Ref<OrbitalSpace> allspace = new OrbitalSpace("", "",
                                                  evecs_ao,
                                                  basis(),
                                                  integral,
                                                  evals,
                                                  0, 0, OrbitalSpace::symmetry);
    std::vector<bool> occmask(nmo);
    for(unsigned int o=0; o<nmo; ++o)
      occmask[o] = (aoccs[o] != 0.0 || boccs[o] != 0.0) ? true : false;
    Ref<OrbitalSpace> occspace = new MaskedOrbitalSpace("", "", allspace, occmask);
    evecs_ao = occspace->coefs();
    evals = occspace->evals();
    nmo = evals.n();

    MOIndexMap o2f = (*allspace << *occspace);
    std::vector<double> aoccs_o(nmo, 1.0);
    std::vector<double> boccs_o(nmo, 1.0);
    for(unsigned int o=0; o<nmo; ++o) {
      const unsigned int oo = o2f[o];
      aoccs_o[o] = aoccs[oo];
      boccs_o[o] = boccs[oo];
    }
    aoccs = aoccs_o;
    boccs = boccs_o;
  }

  // compute active orbital mask
  nmo = evals.n();
  typedef MolecularOrbitalMask<double, RefDiagSCMatrix> FZCMask;
  typedef MolecularOrbitalMask<double, RefDiagSCMatrix, std::greater<double> > FZVMask;
  FZCMask fzcmask(nfzc(), evals);
  FZVMask fzvmask(nfzv(), evals);
  std::vector<bool> actmask(nmo, true);
  // add frozen core and frozen virtuals masks
  std::transform(fzcmask.mask().begin(), fzcmask.mask().end(),
                 fzvmask.mask().begin(), actmask.begin(), std::logical_and<bool>());

  Ref<OrbitalSpaceRegistry> oreg = this->world()->tfactory()->orbital_registry();

  if (scf()->spin_polarized() == false) { // closed-shell
    spinspaces_[Alpha] = new PopulatedOrbitalSpace(oreg, AnySpinCase1, bs, integral, evecs_ao,
                                                   aoccs, actmask, evals, moorder,
                                                   vir_space(), fbrun);
    spinspaces_[Beta] = spinspaces_[Alpha];
  }
  else { // spin-restricted open-shell
    spinspaces_[Alpha] = new PopulatedOrbitalSpace(oreg, Alpha, bs, integral, evecs_ao,
                                                   aoccs, actmask, evals, moorder,
                                                   vir_space(), fbrun);
    spinspaces_[Beta] = new PopulatedOrbitalSpace(oreg, Beta, bs, integral, evecs_ao,
                                                  boccs, actmask, evals, moorder,
                                                  vir_space(), fbrun);
  }
}


void
PsiSCF_RefWavefunction::init_spaces_unrestricted()
{
  // omit unoccupied orbitals?
  const bool omit_uocc = vir_space_.nonnull() && (vir_space_->rank() == 0);
  if (omit_uocc)
    throw FeatureNotImplemented("omit_uocc is not implemented for spin-unrestricted references",
                                __FILE__,__LINE__);

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
  Ref<PsiHSOSHF> hsoshf = dynamic_cast<PsiHSOSHF*>(scf().pointer());
  if (uhf.nonnull() || hsoshf.nonnull()) {
    alpha_evecs = scf()->coefs(Alpha);
    beta_evecs = scf()->coefs(Beta);
    alpha_evals = scf()->evals(Alpha);
    beta_evals = scf()->evals(Beta);
  }
  else {
    throw ProgrammingError("spin-specific spaces not available for this reference function", __FILE__, __LINE__);
  }

  Ref<OrbitalSpaceRegistry> oreg = this->world()->tfactory()->orbital_registry();

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
    spinspaces_[Alpha] = new PopulatedOrbitalSpace(oreg, Alpha, bs, integral,
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
    spinspaces_[Beta] = new PopulatedOrbitalSpace(oreg, Beta, bs, integral,
                                                  beta_evecs,
                                                  bocc, actmask, beta_evals, moorder,
                                                  vir_space(), fbrun);
  }
}

///////////////////////////////////////////////////////////////////

ClassDesc
PsiRASCI_R12RefWavefunction_cd(typeid(PsiRASCI_RefWavefunction),
                     "PsiRASCI_R12RefWavefunction",
                     1,               // version
                     "public R12RefWavefunction", // must match parent
                     0,               // change to create<PsiRASCI_RefWavefunction> if this class is DefaultConstructible
                     0,               // change to 0 if this class is not KeyValConstructible
                     create<PsiRASCI_RefWavefunction>  // change to 0 if this class is not StateInConstructible
                     );

PsiRASCI_RefWavefunction::PsiRASCI_RefWavefunction(const Ref<WavefunctionWorld>& world,
                                                   const Ref<PsiRASCI>& wfn,
                                                   bool spin_restricted,
                                                   unsigned int nfzc,
                                                   unsigned int nfzv,
                                                   bool omit_uocc) :
                                                   RefWavefunction(world,
                                                                      wfn->basis(),
                                                                      wfn->integral()),
                                                   wfn_(wfn),
                                                   spin_restricted_(spin_restricted),
                                                   nfzc_(nfzc),
                                                   nfzv_(nfzv),
                                                   omit_uocc_(omit_uocc) {
  // bring spin_restricted in sync with wfn_
  if (wfn_->spin_polarized() == false) spin_restricted_ = true;
  // if omit_uocc is true, nfzv should be 0
  if (omit_uocc) nfzv = 0;
#if 0
  // make sure that FockBuildRuntime uses same densities as the reference wavefunction
  const double eref = wfn_->energy();
  world->fockbuild_runtime()->set_densities(this->ordm(Alpha), this->ordm(Beta));
#endif
}

PsiRASCI_RefWavefunction::PsiRASCI_RefWavefunction(StateIn& si) : RefWavefunction(si) {
  throw "not implemented";
}

PsiRASCI_RefWavefunction::~PsiRASCI_RefWavefunction() {
  throw "not implemented";
}

void
PsiRASCI_RefWavefunction::save_data_state(StateOut& so) {
  throw "not implemented";
}

void
PsiRASCI_RefWavefunction::init_spaces()
{
  const bool moorder = true;   // order orbitals in the order of increasing energies
  const Ref<GaussianBasisSet> bs = wfn()->basis();
  const Ref<Integral>& integral =  wfn()->integral();
  Ref<PetiteList> plist = integral->petite_list();
  Ref<OrbitalSpace> mospace = wfn()->orbs_sb(Alpha);
  RefSCMatrix evecs_ao = mospace->coefs();
  RefDiagSCMatrix evals = mospace->evals();
  const int nmo = evecs_ao.coldim().n();

  // select occupied orbitals:
  // frozen-core is occupied
  // all RAS1+RAS2 orbitals are always occupied
  // RAS3 are occupied if ras3_max > 0
  using std::vector;
  vector<double> occs(nmo, 0.0);
  vector<bool> occmask(nmo, false);
  const vector<unsigned int> mopi = wfn()->reference()->mopi();
  const vector<unsigned int> frzcpi = wfn()->frozen_docc();
  const vector<unsigned int> ras1 = wfn()->ras1();
  const vector<unsigned int> ras2 = wfn()->ras2();
  const vector<unsigned int> ras3 = wfn()->ras3();
  const unsigned int ras3_max = wfn()->ras3_max();
  const int nirrep = wfn()->nirrep();
  unsigned int mo = 0;
  for(int h=0; h<nirrep; ++h) {
    const unsigned int nmo = mopi.at(h);
    const unsigned int nocc = frzcpi[h] +
                              ras1[h] +
                              ras2[h] +
                              (ras3_max > 0 ? ras3[h] : 0);
    for(int i=0; i<nocc; ++i, ++mo) {
      occs[mo] = 1.0;
      occmask[mo] = true;
    }
    mo += nmo - nocc;
  }

  // omit unoccupied orbitals?
  if (omit_uocc()) {
    Ref<OrbitalSpace> allspace = new OrbitalSpace("", "",
                                                  evecs_ao,
                                                  basis(),
                                                  integral,
                                                  evals,
                                                  0, 0, OrbitalSpace::symmetry);
    Ref<OrbitalSpace> occspace = new MaskedOrbitalSpace("", "", allspace, occmask);
    evecs_ao = occspace->coefs();
    evals = occspace->evals();
    occs = std::vector<double>(evals.n(), 1.0);
  }

  // compute active orbital mask
  typedef MolecularOrbitalMask<double, RefDiagSCMatrix> FZCMask;
  typedef MolecularOrbitalMask<double, RefDiagSCMatrix, std::greater<double> > FZVMask;
  FZCMask fzcmask(nfzc(), evals);
  FZVMask fzvmask(nfzv(), evals);
  std::vector<bool> actmask(nmo, true);
  // add frozen core and frozen virtuals masks
  std::transform(fzcmask.mask().begin(), fzcmask.mask().end(),
                 fzvmask.mask().begin(), actmask.begin(), std::logical_and<bool>());

  Ref<OrbitalSpaceRegistry> oreg = this->world()->tfactory()->orbital_registry();

  // alpha and beta orbitals are the same
  spinspaces_[Alpha] = new PopulatedOrbitalSpace(oreg, AnySpinCase1, bs, integral, evecs_ao,
                                                 occs, actmask, evals, moorder);
  spinspaces_[Beta] = spinspaces_[Alpha];
}

RefSymmSCMatrix
PsiRASCI_RefWavefunction::core_hamiltonian_for_basis(const Ref<GaussianBasisSet> &basis,
                                                     const Ref<GaussianBasisSet> &p_basis)
{
  return wfn()->core_hamiltonian_for_basis(basis, p_basis);
}

RefSymmSCMatrix
PsiRASCI_RefWavefunction::ordm_orbs_sb(SpinCase1 spin) const
{
  RefSymmSCMatrix opdm_full = wfn()->mo_density(spin);
  RefSymmSCMatrix result;
  if (omit_uocc() == false)
    result = opdm_full;
  else {  // if uoccs were omitted map the density from wfn()->orbs_sb() to this->orbs_sb()
    result = opdm_full.kit()->symmmatrix(this->orbs_sb(spin)->dim());
    MOIndexMap o2f = (*(wfn()->orbs_sb(spin)) << *(this->orbs_sb(spin)));
    const unsigned int nocc = this->orbs_sb(spin)->rank();
    for(unsigned int r=0; r<nocc; ++r) {
      unsigned int rr = o2f[r];
      for(unsigned int c=0; c<=r; ++c) {
        unsigned int cc = o2f[c];
        result(r,c) = opdm_full(rr,cc);
      }
    }
  }
  return result;
}

RefSymmSCMatrix
PsiRASCI_RefWavefunction::ordm(SpinCase1 spin) const
{
  RefSymmSCMatrix P = this->ordm_orbs_sb(spin);
  RefSCMatrix C_ao = this->orbs_sb(spin)->coefs();
  RefSymmSCMatrix P_ao = P.kit()->symmmatrix(C_ao.rowdim());
  P_ao.assign(0.0);
  P_ao.accumulate_transform(C_ao, P);
  return P_ao;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
