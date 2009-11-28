//
// refinfo.cc
//
// Copyright (C) 2005 Edward Valeev
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

#include <util/class/scexception.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/mbptr12/refinfo.h>
#include <chemistry/qc/scf/hsosscf.h>
#include <chemistry/qc/basis/obintfactory.h>

using namespace sc;

/////////////

namespace {

  // map occupations from full space to the given space
  RefDiagSCMatrix map_occupations(const Ref<OrbitalSpace>& given,
                                  const Ref<OrbitalSpace>& full,
                                  const RefDiagSCMatrix occs_full) {
    RefDiagSCMatrix occs_given = given->coefs().kit()->diagmatrix(given->coefs().coldim());
    std::vector<unsigned int> given2full = (*full << *given);
    const int n = occs_full.dim().n();
    for(int i=0; i<n; ++i) occs_given(i) = occs_full( given2full[i] );
    return occs_given;
  }

  // converts blocked space into nonblocked space with orbitals ordered in the order of
  // increasing/decreasing orbital energies
  Ref<OrbitalSpace> blocked_to_nonblocked_space(const std::string& id,
                                                const std::string& label,
                                                const Ref<OrbitalSpace>& blocked,
                                                bool eorder_increasing) {
    Ref<OrbitalSpace> nonblocked;
    // occupations are irrelevant for this re-sort
    RefDiagSCMatrix occs = blocked->coefs().kit()->diagmatrix(blocked->coefs().coldim());  occs.assign(0.0);
    if (eorder_increasing) {
      typedef EnergyMOOrder<std::less<double> > Order;
      Order pred;
      nonblocked = new OrderedOrbitalSpace< Order >(id, label,
          blocked->basis(), blocked->integral(), blocked->coefs(), blocked->evals(),
          occs, blocked->orbsym(), pred);
    }
    else {
      typedef EnergyMOOrder<std::greater<double> > Order;
      Order pred;
      nonblocked = new OrderedOrbitalSpace< Order >(id, label,
          blocked->basis(), blocked->integral(), blocked->coefs(), blocked->evals(),
          occs, blocked->orbsym(), pred);
    }
    return nonblocked;
  }

  // given density matrix in AO basis compute its eigenvalues, eigenvectors, and inverse eigenvectors
  void compute_natural_orbitals (const RefSymmSCMatrix& P_ao,
                                 const Ref<PetiteList>& plist,
                                 const Ref<OverlapOrthog>& orthog,
                                 RefDiagSCMatrix& evals,
                                 RefSCMatrix& coefs,
                                 RefSCMatrix& coefs_inv) {
    Ref<GaussianBasisSet> basis = plist->basis();

    // transform density from AO to SO basis
    RefSCMatrix sotoao = plist->sotoao();  // use SO->AO matrix since density is transformed as function, not as operator
    RefSymmSCMatrix P_so = basis->so_matrixkit()->symmmatrix(sotoao.rowdim());
    P_so.assign(0.0);
    P_so.accumulate_transform(sotoao, P_ao);    // P_so = U P_ao U^T

    // transform density from SO to orthogonal SO basis
    RefSCMatrix so2oso = orthog->basis_to_orthog_basis_inverse();
    RefSymmSCMatrix P_oso = P_so.kit()->symmmatrix(orthog->orthog_dim());
    P_oso.assign(0.0);
    P_oso.accumulate_transform(so2oso, P_so, SCMatrix::TransposeTransform);  // P_oso = X^T P_so X

    // diagonalize to obtain natural orbitals (NOs) and occupations
    RefDiagSCMatrix P_evals = P_oso.eigvals();
    RefSCMatrix P_evecs = P_oso.eigvecs();
    P_oso = 0;
    //P_evals.print("density eigenvalues");

    // convert NOs from OSO to AO basis (i.e. compute AO->NO basis and its inverse (NO->AO))
    RefSCMatrix coefs_no =  plist->aotoso() * orthog->basis_to_orthog_basis().t() * P_evecs;
    RefSCMatrix coefs_no_inv = P_evecs.t() * so2oso.t() * sotoao;
    P_evecs = 0;
    so2oso = 0;
    sotoao = 0;

    evals = P_evals;
    coefs = coefs_no;
    coefs_inv = coefs_no_inv;
  }
}

/////////////

const double
PopulatedOrbitalSpace::zero_occupation = 1e-6;

static ClassDesc PopulatedOrbitalSpace_cd(
  typeid(PopulatedOrbitalSpace),"PopulatedOrbitalSpace",1,"virtual public SavableState",
  0, 0, create<PopulatedOrbitalSpace>);


PopulatedOrbitalSpace::PopulatedOrbitalSpace(SpinCase1 spin, const Ref<GaussianBasisSet>& bs,
                                             const Ref<Integral>& integral,
                                             const RefSCMatrix& coefs,
                                             const std::vector<double>& occs,
                                             const std::vector<bool>& active,
                                             const RefDiagSCMatrix& energies,
                                             bool eorder_increasing)
{
  const int nmo = occs.size();
  std::vector<bool> occ_mask(nmo, false);
  std::vector<bool> occ_act_mask(nmo, false);
  std::vector<bool> uocc_mask(nmo, false);
  std::vector<bool> uocc_act_mask(nmo, false);
  for(int i=0; i<nmo; i++) {
    if (fabs(occs[i]) > PopulatedOrbitalSpace::zero_occupation) {
      occ_mask[i] = true;
      occ_act_mask[i] = true && active[i];
    }
    else {
      uocc_mask[i] = true;
      uocc_act_mask[i] = true && active[i];
    }
  }

  const std::string prefix(to_string(spin));
  using std::ostringstream;
  {
    ostringstream oss;
    oss << prefix << " symmetry-blocked MOs";
    std::string id = ParsedOrbitalSpaceKey::key(std::string("p(sym)"),spin);
    orbs_sb_ = new OrbitalSpace(id, oss.str(), coefs, bs, integral, energies, 0, 0, OrbitalSpace::symmetry);
  }
  {
    ostringstream oss;
    oss << prefix << " energy-ordered MOs";
    std::string id = ParsedOrbitalSpaceKey::key(std::string("p"),spin);
    orbs_ = blocked_to_nonblocked_space(id, oss.str(),
                                        orbs_sb_,
                                        eorder_increasing);
  }
  {
    ostringstream oss;
    oss << prefix << " occupied symmetry-blocked MOs";
    std::string id = ParsedOrbitalSpaceKey::key(std::string("m(sym)"),spin);
    occ_sb_ = new MaskedOrbitalSpace(id, oss.str(), orbs_sb_, occ_mask);
  }
  {
    ostringstream oss;
    oss << prefix << " active occupied symmetry-blocked MOs";
    std::string id = ParsedOrbitalSpaceKey::key(std::string("i(sym)"),spin);
    occ_act_sb_ = new MaskedOrbitalSpace(id, oss.str(), orbs_sb_, occ_act_mask);
  }
  {
    ostringstream oss;
    oss << prefix << " occupied MOs";
    std::string id = ParsedOrbitalSpaceKey::key(std::string("m"),spin);
    occ_ = blocked_to_nonblocked_space(id, oss.str(),
                                       occ_sb_,
                                       eorder_increasing);
  }
  {
    ostringstream oss;
    oss << prefix << " active occupied MOs";
    std::string id = ParsedOrbitalSpaceKey::key(std::string("i"),spin);
    occ_act_ = blocked_to_nonblocked_space(id, oss.str(),
                                           occ_act_sb_,
                                           eorder_increasing);
  }
  {
    ostringstream oss;
    oss << prefix << " unoccupied symmetry-blocked MOs";
    std::string id = ParsedOrbitalSpaceKey::key(std::string("e(sym)"),spin);
    uocc_sb_ = new MaskedOrbitalSpace(id, oss.str(), orbs_sb_, uocc_mask);
  }
  {
    ostringstream oss;
    oss << prefix << " active unoccupied symmetry-blocked MOs";
    std::string id = ParsedOrbitalSpaceKey::key(std::string("a(sym)"),spin);
    uocc_act_sb_ = new MaskedOrbitalSpace(id, oss.str(), orbs_sb_, uocc_act_mask);
  }
  {
    ostringstream oss;
    oss << prefix << " unoccupied MOs";
    std::string id = ParsedOrbitalSpaceKey::key(std::string("e"),spin);
    uocc_ = blocked_to_nonblocked_space(id, oss.str(),
                                        uocc_sb_,
                                        eorder_increasing);
  }
  {
    ostringstream oss;
    oss << prefix << " active unoccupied MOs";
    std::string id = ParsedOrbitalSpaceKey::key(std::string("a"),spin);
    uocc_act_ = blocked_to_nonblocked_space(id, oss.str(),
                                            uocc_act_sb_,
                                            eorder_increasing);
  }
}

PopulatedOrbitalSpace::PopulatedOrbitalSpace(StateIn& si) : SavableState(si) {
  orbs_sb_ << SavableState::restore_state(si);
  orbs_ << SavableState::restore_state(si);
  occ_sb_ << SavableState::restore_state(si);
  occ_act_sb_ << SavableState::restore_state(si);
  occ_ << SavableState::restore_state(si);
  occ_act_ << SavableState::restore_state(si);
  uocc_sb_ << SavableState::restore_state(si);
  uocc_act_sb_ << SavableState::restore_state(si);
  uocc_ << SavableState::restore_state(si);
  uocc_act_ << SavableState::restore_state(si);
}

PopulatedOrbitalSpace::~PopulatedOrbitalSpace() {
}

void
PopulatedOrbitalSpace::save_data_state(StateOut& so) {
  SavableState::save_state(orbs_sb_.pointer(),so);
  SavableState::save_state(orbs_.pointer(),so);
  SavableState::save_state(occ_sb_.pointer(),so);
  SavableState::save_state(occ_act_sb_.pointer(),so);
  SavableState::save_state(occ_.pointer(),so);
  SavableState::save_state(occ_act_.pointer(),so);
  SavableState::save_state(uocc_sb_.pointer(),so);
  SavableState::save_state(uocc_act_sb_.pointer(),so);
  SavableState::save_state(uocc_.pointer(),so);
  SavableState::save_state(uocc_act_.pointer(),so);
}

///////////////////////////////////////////////////////////

static ClassDesc RefInfo_cd(
  typeid(RefInfo),"RefInfo",1,"virtual public SavableState",
  0, 0, 0);

RefInfo::RefInfo(const Ref<GaussianBasisSet>& basis,
                 const Ref<Integral>& integral) :
  basis_(basis), integral_(integral) {
  for(int spin=0; spin<NSpinCases1; ++spin) spinspaces_[spin] = 0;
}

RefInfo::RefInfo(StateIn& si) :
  SavableState(si)
{
  basis_ << SavableState::restore_state(si);
  integral_ = Integral::get_default_integral();
  // is the current default Integral compatible with the original factory used to produce this RefInfo?
  Integral::CartesianOrdering o; int io; si.get(io); o = static_cast<Integral::CartesianOrdering>(io);
  if (o != integral_->cartesian_ordering())
    throw InputError("default Integral is incompatible with the Integral used to produce this object",
                     __FILE__,__LINE__);
  integral_->set_basis(basis_);

  for(int spin=0; spin<NSpinCases1; spin++)
    spinspaces_[spin] << SavableState::restore_state(si);
}

RefInfo::~RefInfo()
{
}

void
RefInfo::save_data_state(StateOut& so)
{
  SavableState::save_state(basis_.pointer(),so);
  so.put(static_cast<int>(integral_->cartesian_ordering()));

  for(int spin=0; spin<NSpinCases1; spin++)
    SavableState::save_state(spinspaces_[spin].pointer(),so);
}

void
RefInfo::init() const
{
  if (spinspaces_[Alpha].null()) {
    RefInfo* this_nonconst = const_cast<RefInfo*>(this);
    this_nonconst->init_spaces();
  }
}

const Ref<OrbitalSpace>&
RefInfo::orbs_sb(SpinCase1 s) const
{
  init();
  return spinspaces_[s]->orbs_sb();
}

const Ref<OrbitalSpace>&
RefInfo::orbs(SpinCase1 s) const
{
  init();
  return spinspaces_[s]->orbs();
}

const Ref<OrbitalSpace>&
RefInfo::occ_sb(SpinCase1 s) const
{
  init();
  return spinspaces_[s]->occ_sb();
}

const Ref<OrbitalSpace>&
RefInfo::occ_act_sb(SpinCase1 s) const
{
  init();
  return spinspaces_[s]->occ_act_sb();
}

const Ref<OrbitalSpace>&
RefInfo::occ(SpinCase1 s) const
{
  init();
  return spinspaces_[s]->occ();
}

const Ref<OrbitalSpace>&
RefInfo::occ_act(SpinCase1 s) const
{
  init();
  return spinspaces_[s]->occ_act();
}

const Ref<OrbitalSpace>&
RefInfo::uocc_sb(SpinCase1 s) const
{
  init();
  return spinspaces_[s]->uocc_sb();
}

const Ref<OrbitalSpace>&
RefInfo::uocc_act_sb(SpinCase1 s) const
{
  init();
  return spinspaces_[s]->uocc_act_sb();
}

const Ref<OrbitalSpace>&
RefInfo::uocc(SpinCase1 s) const
{
  init();
  return spinspaces_[s]->uocc();
}

const Ref<OrbitalSpace>&
RefInfo::uocc_act(SpinCase1 s) const
{
  init();
  return spinspaces_[s]->uocc_act();
}

///////////////////////////////////////////////////////////////////

static ClassDesc SlaterDeterminantRefInfo_cd(
  typeid(SlaterDeterminantRefInfo),"SlaterDeterminantRefInfo",1,"public RefInfo",
  0, 0, create<SlaterDeterminantRefInfo>);

SlaterDeterminantRefInfo::SlaterDeterminantRefInfo(const Ref<OneBodyWavefunction>& obwfn,
                                                   bool spin_restricted,
                                                   unsigned int nfzc,
                                                   unsigned int nfzv) :
                                                   RefInfo(obwfn->basis(), obwfn->integral()),
                                                   obwfn_(obwfn),
                                                   spin_restricted_(spin_restricted),
                                                   nfzc_(nfzc),
                                                   nfzv_(nfzv) {
  // bring spin_restricted in sync with obwfn
  if (obwfn_->spin_polarized() == false) spin_restricted_ = true;
  if (obwfn_->spin_unrestricted() == true) spin_restricted_ = false;
}

SlaterDeterminantRefInfo::SlaterDeterminantRefInfo(StateIn& si) : RefInfo(si) {
  obwfn_ << SavableState::restore_state(si);
  si.get(spin_restricted_);
  si.get(nfzc_);
  si.get(nfzv_);
}

SlaterDeterminantRefInfo::~SlaterDeterminantRefInfo() {
}

void
SlaterDeterminantRefInfo::save_data_state(StateOut& so) {
  SavableState::save_state(obwfn_.pointer(), so);
  so.put(spin_restricted_);
  so.put(nfzc_);
  so.put(nfzv_);
}

void
SlaterDeterminantRefInfo::init_spaces()
{
  if (spin_restricted())
    init_spaces_restricted();
  else
    init_spaces_unrestricted();
}

void
SlaterDeterminantRefInfo::init_spaces_restricted()
{
  const Ref<GaussianBasisSet> bs = obwfn()->basis();
  const RefSCMatrix evecs_so = obwfn()->eigenvectors();
  const RefDiagSCMatrix evals = obwfn()->eigenvalues();
  const Ref<Integral>& integral =  obwfn()->integral();
  Ref<PetiteList> plist = integral->petite_list();
  RefSCMatrix evecs_ao = plist->evecs_to_AO_basis(evecs_so);
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
    aoccs[mo] = obwfn()->alpha_occupation(mo);
  }
  if (obwfn()->spin_polarized() == false) { // closed-shell
    spinspaces_[Alpha] = new PopulatedOrbitalSpace(AnySpinCase1, bs, integral, evecs_ao, aoccs, actmask, evals);
    spinspaces_[Beta] = spinspaces_[Alpha];
  }
  else { // spin-restricted open-shell
    spinspaces_[Alpha] = new PopulatedOrbitalSpace(Alpha, bs, integral, evecs_ao, aoccs, actmask, evals);
    vector<double> boccs(nmo);
    for(int mo=0; mo<nmo; mo++) {
      boccs[mo] = obwfn()->beta_occupation(mo);
    }
    spinspaces_[Beta] = new PopulatedOrbitalSpace(Beta, bs, integral, evecs_ao, boccs, actmask, evals);
  }
}


void
SlaterDeterminantRefInfo::init_spaces_unrestricted()
{
  const Ref<GaussianBasisSet> bs = obwfn()->basis();
  const Ref<Integral>& integral = obwfn()->integral();

  using std::vector;
  vector<double> aocc, bocc;
  const int nmo = obwfn()->alpha_eigenvectors().coldim().n();
  for(int mo=0; mo<nmo; mo++) {
    aocc.push_back(obwfn()->alpha_occupation(mo));
    bocc.push_back(obwfn()->beta_occupation(mo));
  }
  Ref<PetiteList> plist = obwfn()->integral()->petite_list();

  RefSCMatrix alpha_evecs, beta_evecs;
  RefDiagSCMatrix alpha_evals, beta_evals;
  // alpha and beta orbitals are available for UHF
  if (obwfn()->spin_unrestricted()) {
    alpha_evecs = obwfn()->alpha_eigenvectors();
    beta_evecs = obwfn()->beta_eigenvectors();
    alpha_evals = obwfn()->alpha_eigenvalues();
    beta_evals = obwfn()->beta_eigenvalues();
  }
  // use semicanonical orbitals for ROHF
  else {
    Ref<HSOSSCF> hsosscf = dynamic_cast<HSOSSCF*>(ref().pointer());
    if (hsosscf.null())
      throw ProgrammingError("spin-specific spaces not available for this reference function", __FILE__, __LINE__);
    alpha_evecs = hsosscf->alpha_semicanonical_eigenvectors();
    beta_evecs = hsosscf->beta_semicanonical_eigenvectors();
    alpha_evals = hsosscf->alpha_semicanonical_eigenvalues();
    beta_evals = hsosscf->beta_semicanonical_eigenvalues();
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
                                                   plist->evecs_to_AO_basis(alpha_evecs),
                                                   aocc, actmask, alpha_evals);
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
                                                  plist->evecs_to_AO_basis(beta_evecs),
                                                  bocc, actmask, beta_evals);
  }
}


///////////////////////////////////////////////////////////////////

static ClassDesc ORDMRefInfo_cd(
  typeid(ORDMRefInfo),"ORDMRefInfo",1,"public RefInfo",
  0, 0, create<ORDMRefInfo>);

ORDMRefInfo::ORDMRefInfo(const Ref<GaussianBasisSet>& basis,
                         const Ref<Integral>& integral,
                         const RefSymmSCMatrix& alpha_1rdm,
                         const RefSymmSCMatrix& beta_1rdm,
                         bool spin_restricted,
                         unsigned int nfzc,
                         bool omit_virtuals) :
                         RefInfo(basis, integral),
                         spin_restricted_(spin_restricted),
                         nfzc_(nfzc),
                         omit_virtuals_(omit_virtuals)
{
  rdm_[Alpha] = alpha_1rdm;
  rdm_[Beta] = beta_1rdm;

  if (rdm_[Alpha] == rdm_[Beta] && spin_restricted_ == false)
    throw ProgrammingError("identical 1-rdms given for alpha and beta spins but spin_restricted = true",__FILE__,__LINE__);
  if (nfzc_ >= basis->nbasis())
    throw ProgrammingError("nfzc > basis set size",__FILE__,__LINE__);
}

ORDMRefInfo::ORDMRefInfo(StateIn& si) : RefInfo(si) {
  int c = 0;
  detail::FromStateIn<RefSymmSCMatrix>::get(rdm_[Alpha], si, c);
  detail::FromStateIn<RefSymmSCMatrix>::get(rdm_[Beta], si, c);
  si.get(spin_restricted_);
  si.get(nfzc_);
  si.get(omit_virtuals_);
}

ORDMRefInfo::~ORDMRefInfo() {
}

void
ORDMRefInfo::save_data_state(StateOut& so) {
  int c = 0;
  detail::ToStateOut<RefSymmSCMatrix>::put(rdm_[Alpha], so, c);
  detail::ToStateOut<RefSymmSCMatrix>::put(rdm_[Alpha], so, c);
  so.put(spin_restricted_);
  so.put(nfzc_);
  so.put(omit_virtuals_);
}

void
ORDMRefInfo::init_spaces() {
  if (spin_restricted())
    init_spaces_restricted();
  else
    init_spaces_unrestricted();
}

void
ORDMRefInfo::init_spaces_restricted() {

  //////////////
  // compute natural orbitals of the total (spin-free) density matrix
  //////////////
  RefSymmSCMatrix P_ao = rdm_[Alpha] + rdm_[Beta];

  // compute overlap in SO basis
  Ref<PetiteList> plist = integral()->petite_list();
  RefSymmSCMatrix S_so = compute_onebody_matrix<&Integral::overlap>(plist);

  // compute orthogonal SO basis transform
  const int debug = 0;
  Ref<OverlapOrthog> orthog = new OverlapOrthog(OverlapOrthog::default_orthog_method(), S_so,
                                                S_so.kit(), OverlapOrthog::default_lindep_tol(),
                                                debug);

  // compute natural orbitals in AO basis (as well as inverse of the coefficient matrix)
  RefDiagSCMatrix P_evals;
  RefSCMatrix coefs_no, coefs_no_inv;
  compute_natural_orbitals(P_ao, plist, orthog, P_evals, coefs_no, coefs_no_inv);

  // compute active orbital mask
  const int nmo = coefs_no_inv.rowdim().n();
  typedef MolecularOrbitalMask<double, RefDiagSCMatrix, std::greater<double> > FZCMask;
  typedef MolecularOrbitalMask<double, RefDiagSCMatrix> FZVMask;
  FZCMask fzcmask(nfzc(), P_evals);
  std::vector<bool> actmask(nmo, true);
  if (omit_virtuals()) {
    // count how many orbitals are unoccupied
    unsigned int nfzv = 0;
    for(int mo=0; mo<nmo; ++mo)
      if (fabs(P_evals(mo)) < PopulatedOrbitalSpace::zero_occupation)
        ++nfzv;
    FZVMask fzvmask(nfzv, P_evals);
    // add frozen core and frozen virtuals masks
    std::transform(fzcmask.mask().begin(), fzcmask.mask().end(),
                   fzvmask.mask().begin(), actmask.begin(), std::logical_and<bool>());
  }
  else // mask out only the frozen core
    actmask = fzcmask.mask();

  // diagonal elements of spin-specific densities in NO basis are the occupation numbers
  RefSymmSCMatrix Pa_no = rdm_[Alpha].kit()->symmmatrix(coefs_no_inv.rowdim());
  Pa_no.assign(0.0);
  Pa_no.accumulate_transform(coefs_no_inv, rdm_[Alpha]);
  using std::vector;
  vector<double> aoccs(nmo);
  for(int mo=0; mo<nmo; mo++) {
    aoccs[mo] = Pa_no(mo, mo);
  }
  Pa_no = 0;
  RefDiagSCMatrix Pa_diag = rdm_[Alpha].kit()->diagmatrix(coefs_no_inv.rowdim());
  Pa_diag.assign(&(aoccs[0]));
  //Pa_diag.print("Alpha occupation numbers");

  const bool no_order = false; // order NOs in decreasing occupation number
  if (spin_polarized() == false) { // closed-shell
    spinspaces_[Alpha] = new PopulatedOrbitalSpace(AnySpinCase1, basis(), integral(), coefs_no,
                                                   aoccs, actmask, Pa_diag, no_order);
    spinspaces_[Beta] = spinspaces_[Alpha];
  }
  else { // spin-restricted open-shell
    spinspaces_[Alpha] = new PopulatedOrbitalSpace(Alpha, basis(), integral(), coefs_no,
                                                   aoccs, actmask, Pa_diag, no_order);
    RefSymmSCMatrix Pb_no = rdm_[Beta].kit()->symmmatrix(coefs_no_inv.rowdim());
    Pb_no.assign(0.0);
    Pb_no.accumulate_transform(coefs_no_inv, rdm_[Beta]);
    vector<double> boccs(nmo);
    for(int mo=0; mo<nmo; mo++) {
      boccs[mo] = Pb_no(mo, mo);
    }
    Pb_no = 0;
    RefDiagSCMatrix Pb_diag = rdm_[Beta].kit()->diagmatrix(coefs_no_inv.rowdim());
    Pb_diag.assign(&(boccs[0]));
    spinspaces_[Beta] = new PopulatedOrbitalSpace(Beta, basis(), integral(), coefs_no,
                                                  boccs, actmask, Pb_diag, no_order);
  }

#if 0
  // reconstruct densities to verify the logic
  {
    RefSCMatrix coefs = spinspaces_[Alpha]->occ()->coefs();
    RefSymmSCMatrix Pr = coefs.kit()->symmmatrix(coefs.rowdim()); Pr.accumulate_symmetric_product(coefs);
    (rdm_[Alpha] - Pr).print("Alpha P - P(reconstruct): should be zero");
  }
  {
    RefSCMatrix coefs = spinspaces_[Beta]->occ()->coefs();
    RefSymmSCMatrix Pr = coefs.kit()->symmmatrix(coefs.rowdim()); Pr.accumulate_symmetric_product(coefs);
    (rdm_[Beta] - Pr).print("Beta P - P(reconstruct): should be zero");
  }
#endif
}

void
ORDMRefInfo::init_spaces_unrestricted() {

  // compute overlap in SO basis
  Ref<PetiteList> plist = integral()->petite_list();
  RefSymmSCMatrix S_so = compute_onebody_matrix<&Integral::overlap>(plist);

  // compute orthogonal SO basis transform
  const int debug = 0;
  Ref<OverlapOrthog> orthog = new OverlapOrthog(OverlapOrthog::default_orthog_method(), S_so,
                                                S_so.kit(), OverlapOrthog::default_lindep_tol(),
                                                debug);

  for(int s=0; s<NSpinCases1; ++s) {
    const SpinCase1 spin = static_cast<SpinCase1>(s);

    // compute natural orbitals in AO basis (as well as inverse of the coefficient matrix)
    RefSymmSCMatrix P_ao = rdm_[s];
    RefDiagSCMatrix P_evals;
    RefSCMatrix coefs_no, coefs_no_inv;
    compute_natural_orbitals(P_ao, plist, orthog, P_evals, coefs_no, coefs_no_inv);

    // compute active orbital mask
    const int nmo = coefs_no_inv.rowdim().n();
    typedef MolecularOrbitalMask<double, RefDiagSCMatrix, std::greater<double> > FZCMask;
    typedef MolecularOrbitalMask<double, RefDiagSCMatrix> FZVMask;
    FZCMask fzcmask(nfzc(), P_evals);
    std::vector<bool> actmask(nmo, true);
    if (omit_virtuals()) {
      // count how many orbitals are unoccupied
      unsigned int nfzv = 0;
      for(int mo=0; mo<nmo; ++mo)
        if (fabs(P_evals(mo)) < PopulatedOrbitalSpace::zero_occupation)
          ++nfzv;
      FZVMask fzvmask(nfzv, P_evals);
      // add frozen core and frozen virtuals masks
      std::transform(fzcmask.mask().begin(), fzcmask.mask().end(),
                     fzvmask.mask().begin(), actmask.begin(), std::logical_and<bool>());
    }
    else // mask out only the frozen core
      actmask = fzcmask.mask();

    using std::vector;
    vector<double> occs(nmo);
    for(int mo=0; mo<nmo; mo++) {
      occs[mo] = P_evals(mo);
    }

    const bool no_order = false; // order NOs in decreasing occupation number
    spinspaces_[spin] = new PopulatedOrbitalSpace(spin, basis(), integral(), coefs_no,
                                                   occs, actmask, P_evals, no_order);

#if 0
    // reconstruct densities to verify the logic
    {
      RefSCMatrix coefs = spinspaces_[spin]->occ()->coefs();
      RefSymmSCMatrix Pr = coefs.kit()->symmmatrix(coefs.rowdim()); Pr.accumulate_symmetric_product(coefs);
      (P_ao - Pr).print(prepend_spincase(spin," P - P(reconstruct): should be zero").c_str());
    }
#endif
  } // end of spincase1 loop
}
