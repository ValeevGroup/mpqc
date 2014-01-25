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
#include <chemistry/qc/mbptr12/pt2r12.h>
#include <chemistry/qc/psi/psiref.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/wfn/spin.h>

using namespace sc;

///////////////////////////////////////////////////////////////////

static ClassDesc PsiSCF_RefWavefunction_cd(
  typeid(PsiSCF_RefWavefunction),"PsiSCF_RefWavefunction",1,"public RefWavefunction",
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
  // spin_restricted is a recommendation only -> make sure it is realizable
  if (scf_->spin_polarized() == false) spin_restricted_ = true;
  Ref<PsiUHF> uhf; uhf << scf_; if (uhf) spin_restricted_ = false;

  if (nfzv > 0 && vir_space)
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

void
PsiSCF_RefWavefunction::print(std::ostream&o) const {
  using std::endl;
  o << indent << "PsiSCF_RefWavefunction:" << endl;
  o << incindent;
    o << indent << "spin_restricted = " << (spin_restricted_ ? "true" : "false") << endl;
    o << indent << "# frozen core   = " << nfzc_ << endl;
    o << indent << "# frozen virt   = " << nfzv_ << endl;
    if (vir_space_) {
      o << indent << "vir_basis:" << endl;
      vir_space_->basis()->print(o);
      o << endl;
    }
    scf_->print(o);
  o << decindent;
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
  RefSymmSCMatrix P;

  Ref<PsiHSOSHF> hsoshf; hsoshf << scf_;

  if (hsoshf && spin_restricted_ == false) { // HSOSHF + spin_unrestricted => must use semicanonical orbitals
    // compute semicanonical densities
    RefSCMatrix C = hsoshf->coefs_semicanonical(s);
    RefDiagSCMatrix P_mo = C.kit()->diagmatrix(C.coldim()); // density in MO basis
    P_mo.assign(0.0);
    const int nmo = P_mo.n();
    for(int mo=0; mo<nmo; ++mo)
      P_mo.set_element(mo, (s == Alpha) ? hsoshf->alpha_occupation(mo)
                                        : hsoshf->beta_occupation(mo) );
    P = C.kit()->symmmatrix(C.rowdim()); P.assign(0.0);
    P.accumulate_transform(C, P_mo);
  }
  else {
    P = (s == Alpha) ? scf()->alpha_ao_density() : scf()->beta_ao_density();
  }
  return P;
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
  const bool omit_uocc = vir_space_ && vir_space_->rank() == 0;
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
  std::vector<ParticleHoleOrbitalAttributes> actmask_a(nmo, ParticleHoleOrbitalAttributes::None);
  std::vector<ParticleHoleOrbitalAttributes> actmask_b(nmo, ParticleHoleOrbitalAttributes::None);
  for(size_t o=0; o<nmo; ++o) {
    actmask_a[o] = aoccs[o] == 0.0 ? ParticleHoleOrbitalAttributes::Particle : ParticleHoleOrbitalAttributes::Hole;
    actmask_b[o] = boccs[o] == 0.0 ? ParticleHoleOrbitalAttributes::Particle : ParticleHoleOrbitalAttributes::Hole;
  }
  typedef MolecularOrbitalMask<double, RefDiagSCMatrix> FZCMask;
  typedef MolecularOrbitalMask<double, RefDiagSCMatrix, std::greater<double> > FZVMask;
  FZCMask fzcmask(nfzc(), evals);
  FZVMask fzvmask(nfzv(), evals);
  // add frozen core and frozen virtuals masks
  for(size_t o=0; o<nmo; ++o) {
    if (not fzcmask[o] || not fzvmask[o]) {
      actmask_a[o] = ParticleHoleOrbitalAttributes::None;
      actmask_b[o] = ParticleHoleOrbitalAttributes::None;
    }
  }

  Ref<OrbitalSpaceRegistry> oreg = this->world()->tfactory()->orbital_registry();

  if (scf()->spin_polarized() == false) { // closed-shell
    spinspaces_[Alpha] = new PopulatedOrbitalSpace(oreg, AnySpinCase1, bs, integral, evecs_ao,
                                                   aoccs, actmask_a, evals, moorder,
                                                   vir_space(), fbrun);
    spinspaces_[Beta] = spinspaces_[Alpha];
  }
  else { // spin-restricted open-shell
    spinspaces_[Alpha] = new PopulatedOrbitalSpace(oreg, Alpha, bs, integral, evecs_ao,
                                                   aoccs, actmask_a, evals, moorder,
                                                   vir_space(), fbrun);
    spinspaces_[Beta] = new PopulatedOrbitalSpace(oreg, Beta, bs, integral, evecs_ao,
                                                  boccs, actmask_b, evals, moorder,
                                                  vir_space(), fbrun);
  }
}


void
PsiSCF_RefWavefunction::init_spaces_unrestricted()
{
  // omit unoccupied orbitals?
  const bool omit_uocc = vir_space_ && (vir_space_->rank() == 0);
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
  if (uhf) {
    alpha_evecs = scf()->coefs(Alpha);
    beta_evecs = scf()->coefs(Beta);
    alpha_evals = scf()->evals(Alpha);
    beta_evals = scf()->evals(Beta);
  }
  else {
    if (hsoshf == 0) // only know what to do with HSOSSCF
      throw ProgrammingError("spin-specific spaces not available for this reference function", __FILE__, __LINE__);

    alpha_evecs = hsoshf->coefs_semicanonical(Alpha);
    beta_evecs = hsoshf->coefs_semicanonical(Beta);
    alpha_evals = hsoshf->evals_semicanonical(Alpha);
    beta_evals = hsoshf->evals_semicanonical(Beta);
  }

  Ref<OrbitalSpaceRegistry> oreg = this->world()->tfactory()->orbital_registry();

  typedef MolecularOrbitalMask<double, RefDiagSCMatrix> FZCMask;
  typedef MolecularOrbitalMask<double, RefDiagSCMatrix, std::greater<double> > FZVMask;
  { // alpha spin
    // compute active orbital mask
    std::vector<ParticleHoleOrbitalAttributes> actmask(nmo, ParticleHoleOrbitalAttributes::None);
    for(size_t o=0; o<nmo; ++o) {
      actmask[o] = aocc[o] == 0.0 ? ParticleHoleOrbitalAttributes::Particle : ParticleHoleOrbitalAttributes::Hole;
    }
    typedef MolecularOrbitalMask<double, RefDiagSCMatrix> FZCMask;
    typedef MolecularOrbitalMask<double, RefDiagSCMatrix, std::greater<double> > FZVMask;
    FZCMask fzcmask(nfzc(), alpha_evals);
    FZVMask fzvmask(nfzv(), alpha_evals);
    // add frozen core and frozen virtuals masks
    for(size_t o=0; o<nmo; ++o) {
      if (not fzcmask[o] || not fzvmask[o]) {
        actmask[o] = ParticleHoleOrbitalAttributes::None;
      }
    }

    spinspaces_[Alpha] = new PopulatedOrbitalSpace(oreg, Alpha, bs, integral,
                                                   alpha_evecs,
                                                   aocc, actmask, alpha_evals, moorder,
                                                   vir_space(), fbrun);
  }
  { // beta spin
    // compute active orbital mask
    std::vector<ParticleHoleOrbitalAttributes> actmask(nmo, ParticleHoleOrbitalAttributes::None);
    for(size_t o=0; o<nmo; ++o) {
      actmask[o] = bocc[o] == 0.0 ? ParticleHoleOrbitalAttributes::Particle : ParticleHoleOrbitalAttributes::Hole;
    }
    typedef MolecularOrbitalMask<double, RefDiagSCMatrix> FZCMask;
    typedef MolecularOrbitalMask<double, RefDiagSCMatrix, std::greater<double> > FZVMask;
    FZCMask fzcmask(nfzc(), alpha_evals);
    FZVMask fzvmask(nfzv(), alpha_evals);
    // add frozen core and frozen virtuals masks
    for(size_t o=0; o<nmo; ++o) {
      if (not fzcmask[o] || not fzvmask[o]) {
        actmask[o] = ParticleHoleOrbitalAttributes::None;
      }
    }

    spinspaces_[Beta] = new PopulatedOrbitalSpace(oreg, Beta, bs, integral,
                                                  beta_evecs,
                                                  bocc, actmask, beta_evals, moorder,
                                                  vir_space(), fbrun);
  }
}

Ref<DensityFittingInfo>
PsiSCF_RefWavefunction::dfinfo() const {
  return use_world_dfinfo() ? const_cast<DensityFittingInfo*>(world()->tfactory()->df_info()) : 0;
}

///////////////////////////////////////////////////////////////////

ClassDesc
PsiRASCI_RefWavefunction_cd(typeid(PsiRASCI_RefWavefunction),
                     "PsiRASCI_RefWavefunction",
                     1,               // version
                     "public RefWavefunction", // must match parent
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
}

PsiRASCI_RefWavefunction::PsiRASCI_RefWavefunction(StateIn& si) : RefWavefunction(si) {
  throw "not implemented";
}

PsiRASCI_RefWavefunction::~PsiRASCI_RefWavefunction() {
}

void
PsiRASCI_RefWavefunction::print(std::ostream&o) const {
  using std::endl;
  o << indent << "PsiRASCI_RefWavefunction:" << endl;
  o << incindent;
    o << indent << "spin_restricted = " << (spin_restricted_ ? "true" : "false") << endl;
    o << indent << "omit_uocc       = " << (omit_uocc_ ? "true" : "false") << endl;
    o << indent << "# frozen core   = " << nfzc_ << endl;
    o << indent << "# frozen virt   = " << nfzv_ << endl;
    wfn_->print(o);
  o << decindent;
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
  int nmo = evecs_ao.coldim().n();

  // select occupied orbitals:
  // frozen-core is occupied
  // all RAS1+RAS2 orbitals are always occupied
  // RAS3 are occupied if ras3_max > 0
  using std::vector;
  vector<double> occs(nmo, 0.0);
  vector<double> rasscf_occs(nmo, 0.0); //to get rasscf orbitals (excluding external virtual orbs)
  vector<bool> occmask(nmo, false);
  const vector<unsigned int> mopi = wfn()->reference()->mopi();
  const vector<unsigned int> frzcpi = wfn()->frozen_docc();
  const vector<unsigned int> ras1 = wfn()->ras1();
  const vector<unsigned int> ras2 = wfn()->ras2();
  const vector<unsigned int> ras3 = wfn()->ras3();
  const unsigned int ras3_max = wfn()->ras3_max();
  const int nirrep = wfn()->nirrep();
  unsigned int mo = 0;
  for(int h=0; h<nirrep; ++h)
  {
    const unsigned int nmo = mopi.at(h);
    const unsigned int nocc = frzcpi[h] +
                              ras1[h] +
                              ras2[h] +
                              (ras3_max > 0 ? ras3[h] : 0);
    const unsigned int rasscf_nocc = frzcpi[h] +
                                  ras1[h] +
                                  ras2[h];
    for(int i=0; i<nocc; ++i, ++mo)
    {
      occs[mo] = 1.0;
      occmask[mo] = true;
      rasscf_occs[mo] = (i<rasscf_nocc)?1.0:0.0;
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
    nmo = evals.n();
  }

  // compute active orbital mask
  std::vector<ParticleHoleOrbitalAttributes> actmask(nmo, ParticleHoleOrbitalAttributes::None);
  for(size_t o=0; o<nmo; ++o) {
    actmask[o] = occs[o] == 0.0 ? ParticleHoleOrbitalAttributes::Particle : ParticleHoleOrbitalAttributes::Hole;
  }
  typedef MolecularOrbitalMask<double, RefDiagSCMatrix> FZCMask;
  typedef MolecularOrbitalMask<double, RefDiagSCMatrix, std::greater<double> > FZVMask;
  FZCMask fzcmask(nfzc(), evals);
  FZVMask fzvmask(nfzv(), evals);
  // add frozen core and frozen virtuals masks
  for(size_t o=0; o<nmo; ++o) {
    if (not fzcmask[o] || not fzvmask[o]) {
      actmask[o] = ParticleHoleOrbitalAttributes::None;
    }
  }

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
PsiRASCI_RefWavefunction::ordm(SpinCase1 spin) const
{
  RefSymmSCMatrix P = wfn()->mo_density(spin);
  RefSCMatrix C_ao = wfn()->orbs_sb(spin)->coefs();
  RefSymmSCMatrix P_ao = P.kit()->symmmatrix(C_ao.rowdim());
  P_ao.assign(0.0);
  P_ao.accumulate_transform(C_ao, P);
  return P_ao;
}

Ref<DensityFittingInfo>
PsiRASCI_RefWavefunction::dfinfo() const {
  return use_world_dfinfo() ? const_cast<DensityFittingInfo*>(world()->tfactory()->df_info()) : 0;
}

/////////////////////////////////////////////////////////////////////////////

Ref<RefWavefunction>
RefWavefunctionFactory::make(const Ref<WavefunctionWorld> & world,
                                const Ref<Wavefunction> & ref,
                                bool spin_restricted,
                                unsigned int nfzc,
                                unsigned int nfzv,
                                Ref<OrbitalSpace> vir_space)
{
  { // PsiSCF
    Ref<PsiSCF> cast; cast << ref;
    if (cast)
      return new PsiSCF_RefWavefunction(world, cast, spin_restricted, nfzc, nfzv, vir_space);
  }
  { // PsiRASCI
    Ref<PsiRASCI> cast; cast << ref;
    if (cast) {
      if (vir_space && vir_space->rank() != 0)
        throw ProgrammingError("PsiRASCI_RefWavefunction can only be used with default virtual space",
                               __FILE__, __LINE__);
      const bool omit_uocc = vir_space;
      return new PsiRASCI_RefWavefunction(world, cast, spin_restricted, nfzc, nfzv, omit_uocc);
    }
  }
  { // OneBodyWavefunction
    Ref<OneBodyWavefunction> cast; cast << ref;
    if (cast)
      return new SD_RefWavefunction(world, cast, spin_restricted, nfzc, nfzv, vir_space);
  }
  throw FeatureNotImplemented("this reference wavefunction cannot be used for correlated methods",
                              __FILE__, __LINE__);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
