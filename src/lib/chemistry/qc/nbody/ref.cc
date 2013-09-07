//
// ref.cc
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

#include <cassert>
#include <util/misc/scexception.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/nbody/ref.h>
#include <chemistry/qc/scf/hsosscf.h>
#include <chemistry/qc/basis/obintfactory.h>
#include <chemistry/qc/wfn/orbitalspace_utils.h>
//#include <chemistry/qc/mbptr12/pt2r12.h>
#include <math/scmat/local.h>
#include <chemistry/qc/lmp2/pop_local.h>
#include <chemistry/qc/lcao/fockbuilder.h>
#include <numeric>

using namespace sc;

/////////////

namespace sc {

  void add_ao_space(const Ref<GaussianBasisSet>& bs,
                    const Ref<Integral>& ints,
                    const Ref<AOSpaceRegistry>& aoreg,
                    const Ref<OrbitalSpaceRegistry> oreg) {
    Ref<Integral> localints = ints->clone();
    const std::string unique_id = new_unique_key(oreg);
    Ref<OrbitalSpace> mu = new AtomicOrbitalSpace(unique_id, std::string(bs->label()) + "(AO)", bs, localints);
    oreg->add(make_keyspace_pair(mu));
    aoreg->add(mu->basis(),mu);
  }
  void remove_ao_space(const Ref<GaussianBasisSet>& bs,
                       const Ref<AOSpaceRegistry>& aoreg,
                       const Ref<OrbitalSpaceRegistry> oreg) {
    Ref<OrbitalSpace> aospace = aoreg->value(bs);
    aoreg->remove(bs);
    oreg->remove(aospace->id());
  }
}

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

Ref<OrbitalSpace>
sc::compute_canonvir_space(const Ref<FockBuildRuntime>& fb_rtime,
                           const Ref<OrbitalSpace>& vir_space,
                           SpinCase1 spin) {
  // add vir_space to the OrbitalSpaceRegistry, if necessary
  const bool need_to_add_virspace = fb_rtime->orbital_registry()->key_exists(vir_space->id()) == false;
  if (need_to_add_virspace)
    fb_rtime->orbital_registry()->add(make_keyspace_pair(vir_space));

  // compute the Fock matrix in VBS space, convert to symmetric form
  const std::string key = ParsedOneBodyIntKey::key(vir_space->id(),vir_space->id(),
                                                   std::string("F"),
                                                   spin);
  RefSCMatrix F = fb_rtime->get(key);
  RefSymmSCMatrix Fs = F.kit()->symmmatrix(F.rowdim());
  Fs.assign(0.0);
  Fs.accumulate_symmetric_sum(F);
  Fs.scale(0.5);

  Ref<OrbitalSpace> result = new OrbitalSpace("e", "canonical symmetry-blocked VBS",
                                              vir_space->coefs()*Fs.eigvecs(),
                                              vir_space->basis(),
                                              vir_space->integral(),
                                              Fs.eigvals(),
                                              0, 0,
                                              OrbitalSpace::symmetry);

  // clean up the mess
  if (need_to_add_virspace)
    fb_rtime->orbital_registry()->remove(vir_space->id());

  return result;
}

/////////////

static ClassDesc PopulatedOrbitalSpace_cd(
  typeid(PopulatedOrbitalSpace),"PopulatedOrbitalSpace",1,"virtual public SavableState",
  0, 0, create<PopulatedOrbitalSpace>);


PopulatedOrbitalSpace::PopulatedOrbitalSpace(const Ref<OrbitalSpaceRegistry>& oreg,
                                             SpinCase1 spin,
                                             const Ref<GaussianBasisSet>& bs,
                                             const Ref<Integral>& integral,
                                             const RefSCMatrix& coefs,
                                             const std::vector<double>& occs,
                                             const std::vector<ParticleHoleOrbitalAttributes>& active,
                                             const RefDiagSCMatrix& energies,
                                             bool eorder_increasing,
                                             Ref<OrbitalSpace> vbs,
                                             Ref<FockBuildRuntime> fbrun
//                                             ,
//                                             std::vector<double> rasscf_occs
                                             ) :
                                             oreg_(oreg)
{
  //
  // validate input
  // 0. sizes agree
  const int nmo = occs.size();
  assert(occs.size() == coefs.ncol());
  assert(bs->nbasis() == coefs.nrow());
  assert(occs.size() == active.size());
  assert(energies.n() == occs.size());
  // 1. occupancies are bounded by 2/1
  assert(fabs(*max_element(occs.begin(), occs.end()) - 1.0) < zero_occupancy() ||
         fabs(*max_element(occs.begin(), occs.end()) - 2.0) < zero_occupancy() );
  assert(*min_element(occs.begin(), occs.end()) >= 0.0);
  // 2. active:
  //   a. can create holes in only occupied orbitals
  //   b. can create particles in only empty or partially occupied orbitals
  for(size_t o=0; o<active.size(); ++o) {
    if (active[o] == ParticleHoleOrbitalAttributes::Hole)
      assert(occs[o] >= zero_occupancy());
    if (active[o] == ParticleHoleOrbitalAttributes::Particle)
      assert(occs[o] <= zero_occupancy());
  }

  //active tells which orbitals are active; the 'masks' are selectors
  std::vector<bool> occ_mask(nmo, false);
  std::vector<bool> occ_act_mask(nmo, false);
  std::vector<bool> uocc_mask(nmo, false);
  std::vector<bool> uocc_act_mask(nmo, false);
  const bool debug = false;
  for(int i=0; i<nmo; i++) {
    if (fabs(occs[i]) > PopulatedOrbitalSpace::zero_occupancy()) {
      occ_mask[i] = true;
      occ_act_mask[i] = (active[i] == ParticleHoleOrbitalAttributes::Hole ||
                         active[i] == ParticleHoleOrbitalAttributes::Any);
    }
    else {
      uocc_mask[i] = true;
      uocc_act_mask[i] = (active[i] == ParticleHoleOrbitalAttributes::Particle ||
                          active[i] == ParticleHoleOrbitalAttributes::Any);
    }
  }

  const std::string prefix(to_string(spin));
  using std::ostringstream;
  {
    ostringstream oss;
    oss << prefix << " symmetry-blocked MOs";
    std::string id = ParsedOrbitalSpaceKey::key(std::string("p"),spin);
    orbs_sb_ = new OrbitalSpace(id, oss.str(), coefs, bs, integral, energies, 0, 0, OrbitalSpace::symmetry);
  }
  {
    ostringstream oss;
    oss << prefix << " energy-ordered MOs";
    std::string id = ParsedOrbitalSpaceKey::key(std::string("p~"),spin);
    orbs_ = blocked_to_nonblocked_space(id, oss.str(),
                                        orbs_sb_,
                                        eorder_increasing);
  }
  {
    ostringstream oss;
    oss << prefix << " occupied symmetry-blocked MOs";
    std::string id = ParsedOrbitalSpaceKey::key(std::string("m"),spin);
    occ_sb_ = new MaskedOrbitalSpace(id, oss.str(), orbs_sb_, occ_mask);
  }

  {
    ostringstream oss;
    oss << prefix << " occupied MOs";
    std::string id = ParsedOrbitalSpaceKey::key(std::string("m~"),spin);
    occ_ = blocked_to_nonblocked_space(id, oss.str(),
                                       occ_sb_,
                                       eorder_increasing);
  }
  {
     ostringstream oss;
     oss << prefix << " active occupied symmetry-blocked MOs";
     std::string id = ParsedOrbitalSpaceKey::key(std::string("i"),spin);
     occ_act_sb_ = new MaskedOrbitalSpace(id, oss.str(), orbs_sb_, occ_act_mask);
  }
  {
    ostringstream oss;
    oss << prefix << " active occupied MOs";
    std::string id = ParsedOrbitalSpaceKey::key(std::string("i~"),spin);
    occ_act_ = blocked_to_nonblocked_space(id, oss.str(),
                                           occ_act_sb_,
                                           eorder_increasing);
  }

  {
    ostringstream oss;
    oss << prefix << " unoccupied symmetry-blocked MOs";
    std::string id = ParsedOrbitalSpaceKey::key(std::string("e"),spin);
    if (vbs.null())
      uocc_sb_ = new MaskedOrbitalSpace(id, oss.str(), orbs_sb_, uocc_mask);
    else {
      if (vbs->rank() > 0) {
        uocc_sb_ = orthog_comp(occ_sb_, vbs, id, "VBS", OverlapOrthog::default_lindep_tol());
        // canonicalize
        uocc_sb_ = compute_canonvir_space(fbrun, uocc_sb_, spin);
      }
      else // empty vbs
        uocc_sb_ = new EmptyOrbitalSpace(id, oss.str(), vbs->basis(), vbs->integral(), OrbitalSpace::symmetry);
    }
  }
  {
    ostringstream oss;
    oss << prefix << " active unoccupied symmetry-blocked MOs";
    std::string id = ParsedOrbitalSpaceKey::key(std::string("a"),spin);
    if (vbs.null())
      uocc_act_sb_ = new MaskedOrbitalSpace(id, oss.str(), orbs_sb_, uocc_act_mask);
    else {
      // make a local uocc_act_mask
      std::vector<bool> _uocc_act_mask(uocc_sb_->rank(),true);
      uocc_act_sb_ = new MaskedOrbitalSpace(id, oss.str(), uocc_sb_, _uocc_act_mask);
    }
  }
  {
    ostringstream oss;
    oss << prefix << " unoccupied MOs";
    std::string id = ParsedOrbitalSpaceKey::key(std::string("e~"),spin);
    uocc_ = blocked_to_nonblocked_space(id, oss.str(),
                                        uocc_sb_,
                                        eorder_increasing);
  }
  {
    ostringstream oss;
    oss << prefix << " active unoccupied MOs";
    std::string id = ParsedOrbitalSpaceKey::key(std::string("a~"),spin);
    uocc_act_ = blocked_to_nonblocked_space(id, oss.str(),
                                            uocc_act_sb_,
                                            eorder_increasing);
  }

  // register all spaces
  Ref<OrbitalSpaceRegistry> idxreg = oreg_;
  idxreg->add(make_keyspace_pair(orbs_sb_));
  idxreg->add(make_keyspace_pair(orbs_));
  idxreg->add(make_keyspace_pair(occ_sb_));
  idxreg->add(make_keyspace_pair(occ_));
  idxreg->add(make_keyspace_pair(uocc_sb_));
  idxreg->add(make_keyspace_pair(uocc_));
  idxreg->add(make_keyspace_pair(occ_act_sb_));
  idxreg->add(make_keyspace_pair(occ_act_));
  idxreg->add(make_keyspace_pair(uocc_act_sb_));
  idxreg->add(make_keyspace_pair(uocc_act_));

#if 0
  orbs_sb_->print_detail();
  occ_sb_->print_detail();
  uocc_sb_->print_detail();
#endif
}



PopulatedOrbitalSpace::PopulatedOrbitalSpace(StateIn& si) : SavableState(si) {
  oreg_ = OrbitalSpaceRegistry::restore_instance(si);
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
  purge();
}

void PopulatedOrbitalSpace::purge() {
  oreg_->remove(orbs_sb_->id());
  oreg_->remove(orbs_->id());
  oreg_->remove(occ_sb_->id());
  oreg_->remove(occ_->id());
  oreg_->remove(occ_act_sb_->id());
  oreg_->remove(occ_act_->id());
  oreg_->remove(uocc_sb_->id());
  oreg_->remove(uocc_->id());
  oreg_->remove(uocc_act_sb_->id());
  oreg_->remove(uocc_act_->id());
}

void
PopulatedOrbitalSpace::save_data_state(StateOut& so) {
  OrbitalSpaceRegistry::save_instance(oreg_, so);
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

static ClassDesc RefWavefunction_cd(
  typeid(RefWavefunction),"RefWavefunction",1,"virtual public SavableState",
  0, 0, 0);

RefWavefunction::RefWavefunction(const Ref<KeyVal>& kv) :
    force_average_AB_rdm1_(false)
{
  world_ << require_dynamic_cast<WavefunctionWorld*>(
      kv->describedclassvalue("world").pointer(),
      "RefWavefunction::RefWavefunction\n"
      );
  if (world_.null())
    throw InputError("RefWavefunction::RefWavefunction -- world is missing");

  basis_ << kv->describedclassvalue("basis").pointer();
  if (basis_.null()) basis_ = world_->wfn()->basis();

  integral_ << kv->describedclassvalue("integral").pointer();
  if (integral_.null()) integral_ = world_->integral()->clone();
  integral_->set_basis(basis_, basis_, basis_, basis_);

  Ref<OneBodyWavefunction> valence_obwfn;
  valence_obwfn << kv->describedclassvalue("valence_obwfn").pointer();
  if (valence_obwfn) {
    RefSCMatrix Cv_so = valence_obwfn->eigenvectors(); // these are in SO basis!
    Ref<GaussianBasisSet> bs = valence_obwfn->basis();
    Ref<Integral> integral = valence_obwfn->integral();
    Ref<PetiteList> plist = integral->petite_list();
    RefSCMatrix Cv_ao = plist->evecs_to_AO_basis(Cv_so);

    valence_orbs_ = new OrbitalSpace("valorbs","valence orbital space",Cv_ao,bs,integral);
  }

  use_world_dfinfo_ = kv->booleanvalue("use_world_df", KeyValValueboolean(false));

  for(int spin=0; spin<NSpinCases1; ++spin) spinspaces_[spin] = 0;
}

RefWavefunction::RefWavefunction(const Ref<WavefunctionWorld>& world,
                                 const Ref<GaussianBasisSet>& basis,
                                 const Ref<Integral>& integral,
                                 bool use_world_df) :
  world_(world), basis_(basis),  integral_(integral->clone()),
  use_world_dfinfo_(use_world_df), force_average_AB_rdm1_(false)
{
  for(int spin=0; spin<NSpinCases1; ++spin) spinspaces_[spin] = 0;
  integral_->set_basis(basis, basis, basis, basis);
}

RefWavefunction::RefWavefunction(StateIn& si) :
  SavableState(si)
{
  world_ << SavableState::restore_state(si);
  basis_ << SavableState::restore_state(si);
  integral_ = Integral::get_default_integral();
  // is the current default Integral compatible with the original factory used to produce this RefWavefunction?
  Integral::CartesianOrdering o; int io; si.get(io); o = static_cast<Integral::CartesianOrdering>(io);
  if (o != integral_->cartesian_ordering())
    throw InputError("default Integral is incompatible with the Integral used to produce this object",
                     __FILE__,__LINE__);
  integral_->set_basis(basis_);
  si.get(use_world_dfinfo_);
  si.get(force_average_AB_rdm1_);

  for(int spin=0; spin<NSpinCases1; spin++)
    spinspaces_[spin] << SavableState::restore_state(si);
}

RefWavefunction::~RefWavefunction()
{
}

void
RefWavefunction::save_data_state(StateOut& so)
{
  SavableState::save_state(world_.pointer(),so);
  SavableState::save_state(basis_.pointer(),so);
  so.put(static_cast<int>(integral_->cartesian_ordering()));
  so.put(use_world_dfinfo_);
  so.put(force_average_AB_rdm1_);

  for(int spin=0; spin<NSpinCases1; spin++)
    SavableState::save_state(spinspaces_[spin].pointer(),so);
}

void
RefWavefunction::set_desired_value_accuracy(double eps) {
  if (eps < this->desired_value_accuracy()) {
    this->reset();
    this->_set_desired_value_accuracy(eps);
  }
}

void
RefWavefunction::init() const
{
  if (spinspaces_[Alpha].null()) {
    RefWavefunction* this_nonconst = const_cast<RefWavefunction*>(this);
    // make sure it's computed first
    const double e = this_nonconst->energy();

    // make sure that FockBuildRuntime uses same density fitting info as this reference
    world_->fockbuild_runtime()->dfinfo(this->dfinfo());

    this_nonconst->init_spaces();//should pay great attention to this!

    // make sure that FockBuildRuntime uses same densities as the reference wavefunction

    RefSymmSCMatrix R = ordm(Alpha);
    if(force_average_AB_rdm1_ == false) // the densities are in AO basis
        world_->fockbuild_runtime()->set_densities(this->ordm(Alpha), this->ordm(Beta));//here computes ordm
    else {
        RefSymmSCMatrix av_rdm = this->ordm(Alpha).copy();
        av_rdm.accumulate(this->ordm(Beta));
        av_rdm.scale(0.5);
        world_->fockbuild_runtime()->set_densities(av_rdm, av_rdm);
    }
  }
}


void
RefWavefunction::obsolete() {
  reset();
}

void
RefWavefunction::purge() {
  obsolete();
}

void
RefWavefunction::reset()
{
  spinspaces_[Alpha] = 0;
  spinspaces_[Beta] = 0;
}

bool
RefWavefunction::desired_value_accuracy_set_to_default() const {
  return true;
}

const Ref<OrbitalSpace>&
RefWavefunction::valence_orbs() const {
  return valence_orbs_;
}

RefSymmSCMatrix
RefWavefunction::ordm_orbs_sb(SpinCase1 spin) const {
  // need to transform density from AO basis to orbs basis; assuming AO density can be obtained from ordm()
  // P' = C^t S P S C
  RefSCMatrix C = orbs_sb(spin)->coefs();
  RefSymmSCMatrix P_ao = this->ordm(spin);
  Ref<PetiteList> plist = integral_->petite_list();
  RefSymmSCMatrix S_so = compute_onebody_matrix<&Integral::overlap>(plist);
  RefSymmSCMatrix S_ao = plist->to_AO_basis(S_so);
  S_so = 0;
  RefSymmSCMatrix SPS_ao = P_ao.kit()->symmmatrix(S_ao.dim()); SPS_ao.assign(0.0);
  SPS_ao.accumulate_transform(S_ao, P_ao);
#if 0
  if (! S_ao.dim()->equiv(P_ao.dim())) { // may need to change the dimension of P_ao to match that of S
    RefSymmSCMatrix P_ao_redim = P_ao.kit()->symmmatrix(S_ao.dim());
    P_ao_redim->convert(P_ao);
    P_ao = P_ao_redim;
    P_ao.print(prepend_spincase(S, "AO density matrix (after redimensioning)").c_str());
  }
#endif
  RefSymmSCMatrix P_mo = C.kit()->symmmatrix(C.coldim());
  P_mo.assign(0.0);
  P_mo.accumulate_transform(C, SPS_ao, SCMatrix::TransposeTransform);
  return P_mo;
}

RefSymmSCMatrix
RefWavefunction::ordm_occ_sb(SpinCase1 spin) const {
  // need to transform density from AO basis to orbs basis; assuming AO density can be obtained from ordm()
  // P' = C^t S P S C
  RefSCMatrix C = occ_sb(spin)->coefs();
  RefSymmSCMatrix P_ao = this->ordm(spin);
  Ref<PetiteList> plist = integral_->petite_list();
  RefSymmSCMatrix S_so = compute_onebody_matrix<&Integral::overlap>(plist);
  RefSymmSCMatrix S_ao = plist->to_AO_basis(S_so);
  S_so = 0;
  RefSymmSCMatrix SPS_ao = P_ao.kit()->symmmatrix(S_ao.dim()); SPS_ao.assign(0.0);
  SPS_ao.accumulate_transform(S_ao, P_ao);
#if 0
  if (! S_ao.dim()->equiv(P_ao.dim())) { // may need to change the dimension of P_ao to match that of S
    RefSymmSCMatrix P_ao_redim = P_ao.kit()->symmmatrix(S_ao.dim());
    P_ao_redim->convert(P_ao);
    P_ao = P_ao_redim;
    P_ao.print(prepend_spincase(S, "AO density matrix (after redimensioning)").c_str());
  }
#endif
  RefSymmSCMatrix P_mo = C.kit()->symmmatrix(C.coldim());
  P_mo.assign(0.0);
  P_mo.accumulate_transform(C, SPS_ao, SCMatrix::TransposeTransform);
  return P_mo;
}

namespace {
  SpinCase1 valid_spincase(SpinCase1 s) {
    return s == AnySpinCase1 ? Alpha : s;
  }
}

const Ref<OrbitalSpace>&
RefWavefunction::orbs_sb(SpinCase1 s) const
{
  init();
  s = valid_spincase(s);
  return spinspaces_[s]->orbs_sb();
}

const Ref<OrbitalSpace>&
RefWavefunction::occ_sb(SpinCase1 s) const
{
  init();
  s = valid_spincase(s);
  return spinspaces_[s]->occ_sb();
}

const Ref<OrbitalSpace>&
RefWavefunction::occ_act_sb(SpinCase1 s) const
{
  init();
  s = valid_spincase(s);
  return spinspaces_[s]->occ_act_sb();
}

const Ref<OrbitalSpace>&
RefWavefunction::uocc_sb(SpinCase1 s) const
{
  init();
  s = valid_spincase(s);
  return spinspaces_[s]->uocc_sb();
}

const Ref<OrbitalSpace>&
RefWavefunction::uocc_act_sb(SpinCase1 s) const
{
  init();
  s = valid_spincase(s);
  return spinspaces_[s]->uocc_act_sb();
}

const Ref<OrbitalSpace>&
RefWavefunction::orbs(SpinCase1 s) const
{
  init();
  s = valid_spincase(s);
  return spinspaces_[s]->orbs();
}
const Ref<OrbitalSpace>&
RefWavefunction::occ(SpinCase1 s) const
{
  init();
  s = valid_spincase(s);
  return spinspaces_[s]->occ();
}

const Ref<OrbitalSpace>&
RefWavefunction::occ_act(SpinCase1 s) const
{
  init();
  s = valid_spincase(s);
  return spinspaces_[s]->occ_act();
}



const Ref<OrbitalSpace>&
RefWavefunction::uocc(SpinCase1 s) const
{
  init();
  s = valid_spincase(s);
  return spinspaces_[s]->uocc();
}

const Ref<OrbitalSpace>&
RefWavefunction::uocc_act(SpinCase1 s) const
{
  init();
  s = valid_spincase(s);
  return spinspaces_[s]->uocc_act();
}

///////////////////////////////////////////////////////////////////

static ClassDesc SD_RefWavefunction_cd(
  typeid(SD_RefWavefunction),"SD_RefWavefunction",1,"public RefWavefunction",
  0, create<SD_RefWavefunction>, create<SD_RefWavefunction>);

SD_RefWavefunction::SD_RefWavefunction(const Ref<KeyVal>& keyval) :
    RefWavefunction(keyval)
{
  obwfn_ << require_dynamic_cast<OneBodyWavefunction*>(
      keyval->describedclassvalue("obwfn").pointer(),
      "SD_RefWavefunction::SD_RefWavefunction\n"
  );
  if (obwfn_.null())
    throw InputError("SD_RefWavefunction::SD_RefWavefunction -- obwfn is missing");

  // spin_restricted is a recommendation only -> make sure it is realizable
  bool default_spin_restricted;
  if (obwfn_->spin_polarized() == false) default_spin_restricted = true;
  if (obwfn_->spin_unrestricted() == true) default_spin_restricted = false;
  spin_restricted_ = keyval->booleanvalue("spin_restricted", KeyValValueboolean(default_spin_restricted));

  nfzc_ = keyval->sizevalue("nfzc", KeyValValuesize(0));
  nfzv_ = keyval->sizevalue("nfzv", KeyValValuesize(0));

  Ref<GaussianBasisSet> vbs; vbs << keyval->describedclassvalue("vir_basis").pointer();
  if (vbs.null())
    vir_space_ = 0;
  else { // check if vbs is known already
    Ref<OrbitalSpaceRegistry> oreg = this->world()->moints_runtime4()->factory()->orbital_registry();
    Ref<AOSpaceRegistry> aoreg = this->world()->moints_runtime4()->factory()->ao_registry();
    if (not aoreg->key_exists(vbs)) { // add AO space, if not known
      add_ao_space(vbs, this->integral(), aoreg, oreg);
    }
    // vbs is labeled e(sym)
    if (oreg->key_exists("e(sym)")) {
      vir_space_ = oreg->value("e(sym)");
      if (not vir_space_->basis()->equiv(vbs))
        throw ProgrammingError("SD_RefWavefunction::SD_RefWavefunction -- vir_space already known and does not match vir_basis",
                               __FILE__, __LINE__);
    }
    else {
      int nlindep_vir = -1; // will compute the number of linear dependencies automatically
      vir_space_ = orthogonalize("e(sym)", "VBS", vbs,
                                 this->integral(),
                                 obwfn_->orthog_method(), obwfn_->lindep_tol(),
                                 nlindep_vir);
    }
  }

  KeyValValuestring def_occ_orbitals("canonical");
  occ_orbitals_ = keyval->stringvalue("occ_orbitals", def_occ_orbitals);
  if (occ_orbitals_ != "pipek-mezey"
      && occ_orbitals_ != "canonical") {
      throw InputError("invalid keyword value",
                       __FILE__, __LINE__, "occ_orbitals", occ_orbitals_.c_str(),
                       this->class_desc());
    }
  if (occ_orbitals_ != "canonical" && this->obwfn()->molecule()->point_group()->order() != 1)
    throw InputError("localized orbitals can only be used in C1 symmetry",
                     __FILE__, __LINE__, "occ_orbitals", occ_orbitals_.c_str(), this->class_desc());
  if (occ_orbitals_ != "canonical" && this->obwfn()->spin_polarized() == true)
    throw InputError("localized occupied orbitals can only be used for closed-shell molecules",
                     __FILE__, __LINE__, "occ_orbitals", occ_orbitals_.c_str(), this->class_desc());

  if (nfzv_ > 0 && vir_space_.nonnull()) {
    std::ostringstream oss;
    oss << nfzv_;
    throw InputError("when VBS is given nfzv must be 0",
                     __FILE__,__LINE__,
                     "nfzv", oss.str().c_str());
  }
}

SD_RefWavefunction::SD_RefWavefunction(const Ref<WavefunctionWorld>& world,
                                             const Ref<OneBodyWavefunction>& obwfn,
                                             bool spin_restricted,
                                             unsigned int nfzc,
                                             unsigned int nfzv,
                                             Ref<OrbitalSpace> vir_space,
                                             std::string occ_orbitals) :
                                             RefWavefunction(world,
                                                                obwfn->basis(),
                                                                obwfn->integral()),
                                             obwfn_(obwfn),
                                             spin_restricted_(spin_restricted),
                                             nfzc_(nfzc),
                                             nfzv_(nfzv),
                                             vir_space_(vir_space),
                                             occ_orbitals_(occ_orbitals)
{
  // spin_restricted is a recommendation only -> make sure it is realizable
  if (obwfn_->spin_polarized() == false) spin_restricted_ = true;
  if (obwfn_->spin_unrestricted() == true) spin_restricted_ = false;

  if (nfzv > 0 && vir_space.nonnull())
    throw ProgrammingError("when VBS is given nfzv must be 0",__FILE__,__LINE__);
}

SD_RefWavefunction::SD_RefWavefunction(StateIn& si) : RefWavefunction(si) {
  obwfn_ << SavableState::restore_state(si);
  vir_space_ << SavableState::restore_state(si);
  si.get(spin_restricted_);
  si.get(nfzc_);
  si.get(nfzv_);
  si.get(occ_orbitals_);
}

SD_RefWavefunction::~SD_RefWavefunction() {
}

void
SD_RefWavefunction::save_data_state(StateOut& so) {
  RefWavefunction::save_data_state(so);
  SavableState::save_state(obwfn_.pointer(), so);
  SavableState::save_state(vir_space_.pointer(), so);
  so.put(spin_restricted_);
  so.put(nfzc_);
  so.put(nfzv_);
  so.put(occ_orbitals_);
}

void
SD_RefWavefunction::print(std::ostream&o) const {
  using std::endl;
  o << indent << "SD_RefWavefunction:" << endl;
  o << incindent;
    o << indent << "spin_restricted = " << (spin_restricted_ ? "true" : "false") << endl;
    o << indent << "# frozen core   = " << nfzc_ << endl;
    o << indent << "# frozen virt   = " << nfzv_ << endl;
    o << indent << "occ_orbitals    = " << occ_orbitals_ << endl;
    if (vir_space_.nonnull()) {
      o << indent << "vir_basis:" << endl;
      vir_space_->basis()->print(o);
      o << endl;
    }
  o << decindent;
}

void
SD_RefWavefunction::purge() {
  obwfn_->purge();
  RefWavefunction::purge();
}

RefSymmSCMatrix
SD_RefWavefunction::ordm(SpinCase1 s) const {
  s = valid_spincase(s);
  if (spin_polarized() == false) s = Alpha;
  if (s == Alpha) return obwfn()->alpha_ao_density();
  return obwfn()->beta_ao_density();
}

RefSymmSCMatrix
SD_RefWavefunction::core_hamiltonian_for_basis(const Ref<GaussianBasisSet> &basis,
                                               const Ref<GaussianBasisSet> &p_basis) {
  return obwfn()->core_hamiltonian_for_basis(basis, p_basis);
}

void
SD_RefWavefunction::init_spaces()
{
  if (spin_restricted())
    init_spaces_restricted();
  else
    init_spaces_unrestricted();
}

void
SD_RefWavefunction::init_spaces_restricted()
{
  const bool moorder = true;   // order orbitals in the order of increasing energies
  Ref<FockBuildRuntime> fbrun = this->world()->fockbuild_runtime();
  RefSymmSCMatrix P_ao = obwfn()->ao_density(); P_ao.scale(0.5);
  fbrun->set_densities(P_ao, P_ao);
  const Ref<GaussianBasisSet> bs = obwfn()->basis();
  RefSCMatrix evecs_so = obwfn()->eigenvectors();
  RefDiagSCMatrix evals = obwfn()->eigenvalues();
  const Ref<Integral>& integral =  obwfn()->integral();
  Ref<PetiteList> plist = integral->petite_list();
  RefSCMatrix evecs_ao = plist->evecs_to_AO_basis(evecs_so);
  int nmo = evecs_ao.coldim().n();

  // optionally, localize the occupied orbitals
  if (occ_orbitals_ == "pipek-mezey") {
    RefSymmSCMatrix S_ao = detail::overlap(bs, integral);
    PipekMezeyLocalization localizer(obwfn(), nfzc_, S_ao);
    //evecs_ao.print("Orbitals before localization");
    RefSCMatrix evecs_occ_ao = localizer.compute_orbitals();
    const int nocc_act = evecs_occ_ao.ncol();
    const int nocc = nocc_act + nfzc_;
    //evecs_occ_ao.print("Active occupied orbitals after localization");
    evecs_ao.assign_subblock(evecs_occ_ao,
                             0, evecs_occ_ao.nrow()-1,
                             nfzc_, nocc - 1);
    //evecs_ao.print("Orbitals after localization");

    // let's test them now!
#if 0
    { // metric
      RefSymmSCMatrix S_mo = evecs_ao.kit()->symmmatrix(evecs_ao.coldim());
      S_mo.assign(0.0);
      S_mo.accumulate_transform(evecs_ao, S_ao, SCMatrix::TransposeTransform);
      S_mo.print("MO overlap");
    }
#endif
    { // recompute the Fock matrix and pick the diagonal elements as the eigenvalues for the localized orbitals
      Ref<FockBuildRuntime> fbrun = this->world()->fockbuild_runtime();
      Ref<OrbitalSpace> aospace = this->world()->tfactory()->orbital_registry()->value("mu");
      std::string F_ao_key = ParsedOneBodyIntKey::key(aospace->id(), aospace->id(), "F", AnySpinCase1);
      RefSCMatrix F_ao_rect = fbrun->get(F_ao_key);
      RefSCMatrix F_mo = evecs_ao.t() * F_ao_rect * evecs_ao;
      for(int o=nfzc_; o<nocc; ++o)
        evals(o) = F_mo.get_element(o,o);
    }
  }

  using std::vector;
  vector<double> aoccs(nmo);
  vector<double> boccs(nmo);
  for(int mo=0; mo<nmo; mo++) {
    aoccs[mo] = obwfn()->alpha_occupation(mo);
    boccs[mo] = obwfn()->beta_occupation(mo);
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
  std::vector<ParticleHoleOrbitalAttributes> actmask_a(nmo, ParticleHoleOrbitalAttributes::None);
  std::vector<ParticleHoleOrbitalAttributes> actmask_b(nmo, ParticleHoleOrbitalAttributes::None);
  nmo = evals.n();
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
  if (obwfn()->spin_polarized() == false) { // closed-shell
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
SD_RefWavefunction::init_spaces_unrestricted()
{
  // omit unoccupied orbitals?
  const bool omit_uocc = vir_space_.nonnull() && (vir_space_->rank() == 0);
  if (omit_uocc)
    throw FeatureNotImplemented("omit_uocc is not implemented for spin-unrestricted references",
                                __FILE__,__LINE__);

  const bool moorder = true;   // order orbitals in the order of increasing energies
  Ref<FockBuildRuntime> fbrun = this->world()->fockbuild_runtime();
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
    Ref<HSOSSCF> hsosscf = dynamic_cast<HSOSSCF*>(obwfn().pointer());
    if (hsosscf.null())
      throw ProgrammingError("spin-specific spaces not available for this reference function", __FILE__, __LINE__);
    alpha_evecs = hsosscf->alpha_semicanonical_eigenvectors();
    beta_evecs = hsosscf->beta_semicanonical_eigenvectors();
    alpha_evals = hsosscf->alpha_semicanonical_eigenvalues();
    beta_evals = hsosscf->beta_semicanonical_eigenvalues();
  }

  Ref<OrbitalSpaceRegistry> oreg = this->world()->tfactory()->orbital_registry();
  typedef MolecularOrbitalMask<double, RefDiagSCMatrix> FZCMask;
  typedef MolecularOrbitalMask<double, RefDiagSCMatrix, std::greater<double> > FZVMask;
  { // alpha spin
    std::vector<ParticleHoleOrbitalAttributes> actmask(nmo, ParticleHoleOrbitalAttributes::None);
    for(size_t o=0; o<nmo; ++o) {
      actmask[o] = aocc[o] == 0.0 ? ParticleHoleOrbitalAttributes::Particle : ParticleHoleOrbitalAttributes::Hole;
    }
    // compute active orbital mask
    FZCMask fzcmask(nfzc(), alpha_evals);
    FZVMask fzvmask(nfzv(), alpha_evals);
    // add frozen core and frozen virtuals masks
    for(size_t o=0; o<nmo; ++o) {
      if (not fzcmask[o] || not fzvmask[o]) {
        actmask[o] = ParticleHoleOrbitalAttributes::None;
      }
    }
    spinspaces_[Alpha] = new PopulatedOrbitalSpace(oreg, Alpha, bs, integral,
                                                   plist->evecs_to_AO_basis(alpha_evecs),
                                                   aocc, actmask, alpha_evals, moorder,
                                                   vir_space(), fbrun);
  }
  { // beta spin
    std::vector<ParticleHoleOrbitalAttributes> actmask(nmo, ParticleHoleOrbitalAttributes::None);
    for(size_t o=0; o<nmo; ++o) {
      actmask[o] = bocc[o] == 0.0 ? ParticleHoleOrbitalAttributes::Particle : ParticleHoleOrbitalAttributes::Hole;
    }
    // compute active orbital mask
    FZCMask fzcmask(nfzc(), beta_evals);
    FZVMask fzvmask(nfzv(), beta_evals);
    // add frozen core and frozen virtuals masks
    for(size_t o=0; o<nmo; ++o) {
      if (not fzcmask[o] || not fzvmask[o]) {
        actmask[o] = ParticleHoleOrbitalAttributes::None;
      }
    }
    spinspaces_[Beta] = new PopulatedOrbitalSpace(oreg, Beta, bs, integral,
                                                  plist->evecs_to_AO_basis(beta_evecs),
                                                  bocc, actmask, beta_evals, moorder,
                                                  vir_space(), fbrun);
  }
}

Ref<DensityFittingInfo>
SD_RefWavefunction::dfinfo() const {
  Ref<DensityFittingInfo> result;
  if (use_world_dfinfo())
    result = const_cast<DensityFittingInfo*>(world()->tfactory()->df_info());
  else {
    Ref<SCF> scf_ptr; scf_ptr << this->obwfn();
    result = (scf_ptr.nonnull()) ? scf_ptr->dfinfo() : 0;
  }
  return result;
}

///////////////////////////////////////////////////////////////////

static ClassDesc Extern_RefWavefunction_cd(
  typeid(Extern_RefWavefunction),"Extern_RefWavefunction",1,"public RefWavefunction",
  0, 0, create<Extern_RefWavefunction>);

Extern_RefWavefunction::Extern_RefWavefunction(const Ref<WavefunctionWorld>& world,
                         const Ref<GaussianBasisSet>& basis,
                         const Ref<Integral>& integral,
                         const RefSCMatrix& orbs,
                         const std::vector<unsigned int>& orbsym,
                         const RefSymmSCMatrix& alpha_1rdm,
                         const RefSymmSCMatrix& beta_1rdm,
                         unsigned int nocc,
                         unsigned int nfzc,
                         unsigned int nfzv,
                         bool omit_uocc) :
                         RefWavefunction(world, basis, integral),
                         nfzc_(nfzc),
                         nfzv_(nfzv),
                         omit_uocc_(omit_uocc)
{
  const unsigned int norbs = orbs.coldim().n();
  assert(nocc >= nfzc);
  assert(norbs - nocc >= nfzv);

  // are densities idempotent?
  ordm_idempotent_ = true;
  {
    RefSymmSCMatrix Pa_2 = alpha_1rdm.kit()->symmmatrix(alpha_1rdm.dim());
    RefSymmSCMatrix I = alpha_1rdm.kit()->symmmatrix(alpha_1rdm.dim());
    Pa_2.assign(0.0);
    I.assign(0.0);
    I->shift_diagonal(1.0);
    Pa_2.accumulate_transform(alpha_1rdm, I); // ?
    Ref<SCElementKNorm> op = new SCElementKNorm;
    (Pa_2 - alpha_1rdm)->element_op(op.pointer());
    if (op->result() > DBL_EPSILON)
      ordm_idempotent_ = false;
    if (beta_1rdm == alpha_1rdm && ordm_idempotent_ == true) {
      RefSymmSCMatrix Pb_2 = beta_1rdm.kit()->symmmatrix(beta_1rdm.dim());
      RefSymmSCMatrix I = beta_1rdm.kit()->symmmatrix(beta_1rdm.dim());
      Pb_2.assign(0.0);
      I.assign(0.0); I->shift_diagonal(1.0);
      Pb_2.accumulate_transform(beta_1rdm, I);
      Ref<SCElementKNorm> op = new SCElementKNorm;
      (Pb_2 - beta_1rdm)->element_op(op.pointer());
      if (op->result() > DBL_EPSILON)
        ordm_idempotent_ = false;
    }
  }

  // how many electrons?
  const double nelectron_alpha = alpha_1rdm.trace();
  const double nelectron_beta = beta_1rdm.trace();
  const double nelectron_double = nelectron_alpha + nelectron_beta;
  assert(nelectron_double > 0.0);
  nelectron_ = round(nelectron_double);
  if (fabs((double)nelectron_ - nelectron_double) > 0.01)
    throw InputError("Extern_RefWavefunction: input density trace is not integer, something is seriously wrong",
                     __FILE__, __LINE__);
  const double magmom_double = nelectron_alpha - nelectron_beta;
  magmom_ = round(magmom_double);
  if (fabs(magmom_ - magmom_double) > 0.01 )
    throw InputError("Extern_RefWavefunction: input density magnetic moment is not half-integer, something is seriously wrong",
                     __FILE__, __LINE__);

  // convert densities to AO basis
  { // Alpha spin
    RefSymmSCMatrix alpha_1rdm_ao = orbs.kit()->symmmatrix(orbs->rowdim());
    alpha_1rdm_ao.assign(0.0);
    alpha_1rdm_ao.accumulate_transform(orbs, alpha_1rdm);
    rdm_[Alpha] = alpha_1rdm_ao;
  }
  if (alpha_1rdm == beta_1rdm)
    rdm_[Beta] = rdm_[Alpha];
  else {
    // Beta spin
    RefSymmSCMatrix beta_1rdm_ao = orbs.kit()->symmmatrix(orbs->rowdim());
    beta_1rdm_ao.assign(0.0);
    beta_1rdm_ao.accumulate_transform(orbs, beta_1rdm);
    rdm_[Beta] = beta_1rdm_ao;
  }

  force_average_AB_rdm1_ = true; // this is meant to always be used with spin-free algorithms
}

Extern_RefWavefunction::Extern_RefWavefunction(const Ref<WavefunctionWorld>& world,
                         const Ref<GaussianBasisSet>& basis,
                         const Ref<Integral>& integral,
                         const RefSCMatrix& orbs,
                         const std::vector<unsigned int>& orbsym,
                         const RefSymmSCMatrix& alpha_1rdm,
                         const RefSymmSCMatrix& beta_1rdm,
                         std::vector<unsigned int> occpi,
                         std::vector<unsigned int> fzcpi,
                         std::vector<unsigned int> fzvpi,
                         std::vector<unsigned int> holepi,
                         std::vector<unsigned int> partpi,
                         bool omit_uocc) :
                         RefWavefunction(world, basis, integral),
                         omit_uocc_(omit_uocc)
{
  if (omit_uocc_)
    throw FeatureNotImplemented("omit_uocc=true is not yet implemented",
                                __FILE__, __LINE__, this->class_desc());

  const unsigned int nfzc = std::accumulate(fzcpi.begin(), fzcpi.end(), 0.0);
  const unsigned int nocc = std::accumulate(occpi.begin(), occpi.end(), 0.0);
  const unsigned int nfzv = std::accumulate(fzvpi.begin(), fzvpi.end(), 0.0);
  nfzc_ = nfzc;
  nfzv_ = nfzv;
  const unsigned int norbs = orbs.coldim().n();
  assert(nocc >= nfzc);
  assert(norbs - nocc >= nfzv);

  // are densities idempotent?
  ordm_idempotent_ = true;
  {
    RefSymmSCMatrix Pa_2 = alpha_1rdm.kit()->symmmatrix(alpha_1rdm.dim());
    RefSymmSCMatrix I = alpha_1rdm.kit()->symmmatrix(alpha_1rdm.dim());
    Pa_2.assign(0.0);
    I.assign(0.0);
    I->shift_diagonal(1.0);
    Pa_2.accumulate_transform(alpha_1rdm, I); // ?
    Ref<SCElementKNorm> op = new SCElementKNorm;
    (Pa_2 - alpha_1rdm)->element_op(op.pointer());
    if (op->result() > DBL_EPSILON)
      ordm_idempotent_ = false;
    if (beta_1rdm == alpha_1rdm && ordm_idempotent_ == true) {
      RefSymmSCMatrix Pb_2 = beta_1rdm.kit()->symmmatrix(beta_1rdm.dim());
      RefSymmSCMatrix I = beta_1rdm.kit()->symmmatrix(beta_1rdm.dim());
      Pb_2.assign(0.0);
      I.assign(0.0); I->shift_diagonal(1.0);
      Pb_2.accumulate_transform(beta_1rdm, I);
      Ref<SCElementKNorm> op = new SCElementKNorm;
      (Pb_2 - beta_1rdm)->element_op(op.pointer());
      if (op->result() > DBL_EPSILON)
        ordm_idempotent_ = false;
    }
  }

  // how many electrons?
  const double nelectron_alpha = alpha_1rdm.trace();
  const double nelectron_beta = beta_1rdm.trace();
  const double nelectron_double = nelectron_alpha + nelectron_beta;
  assert(nelectron_double > 0.0);
  nelectron_ = round(nelectron_double);
  if (fabs((double)nelectron_ - nelectron_double) > 0.01)
    throw InputError("Extern_RefWavefunction: input density trace is not integer, something is seriously wrong",
                     __FILE__, __LINE__);
  const double magmom_double = nelectron_alpha - nelectron_beta;
  magmom_ = round(magmom_double);
  if (fabs(magmom_ - magmom_double) > 0.01 )
    throw InputError("Extern_RefWavefunction: input density magnetic moment is not half-integer, something is seriously wrong",
                     __FILE__, __LINE__);

  // convert densities to AO basis
  { // Alpha spin
    RefSymmSCMatrix alpha_1rdm_ao = orbs.kit()->symmmatrix(orbs->rowdim());
    alpha_1rdm_ao.assign(0.0);
    alpha_1rdm_ao.accumulate_transform(orbs, alpha_1rdm);
    rdm_[Alpha] = alpha_1rdm_ao;
  }
  if (alpha_1rdm == beta_1rdm)
    rdm_[Beta] = rdm_[Alpha];
  else {
    // Beta spin
    RefSymmSCMatrix beta_1rdm_ao = orbs.kit()->symmmatrix(orbs->rowdim());
    beta_1rdm_ao.assign(0.0);
    beta_1rdm_ao.accumulate_transform(orbs, beta_1rdm);
    rdm_[Beta] = beta_1rdm_ao;
  }

 // // this object will never become obsolete, so can compute it right now
  force_average_AB_rdm1_ = true; // this is meant to always be used with spin-free algorithms
  init(orbs, orbsym, occpi, fzcpi, fzvpi, holepi, partpi);
}


Extern_RefWavefunction::Extern_RefWavefunction(StateIn& si) : RefWavefunction(si) {
  int c = 0;
  FromStateIn(rdm_[Alpha], si, c);
  FromStateIn(rdm_[Beta], si, c);
  si.get(nfzc_);
  si.get(nfzv_);
  si.get(nelectron_);
  si.get(magmom_);
  si.get(omit_uocc_);
  si.get(ordm_idempotent_);
}

Extern_RefWavefunction::~Extern_RefWavefunction() {
}

void
Extern_RefWavefunction::save_data_state(StateOut& so) {
  RefWavefunction::save_data_state(so);
  int c = 0;
  ToStateOut(rdm_[Alpha], so, c);
  ToStateOut(rdm_[Alpha], so, c);
  so.put(nfzc_);
  so.put(nfzv_);
  so.put(nelectron_);
  so.put(magmom_);
  so.put(omit_uocc_);
  so.put(ordm_idempotent_);
}

void
Extern_RefWavefunction::print(std::ostream&o) const {
  using std::endl;
  o << indent << "Extern_RefWavefunction:" << endl;
  o << incindent;
    o << indent << "omit_uocc       = " << (omit_uocc_ ? "true" : "false") << endl;
    o << indent << "# frozen core   = " << nfzc_ << endl;
    o << indent << "# frozen virt   = " << nfzv_ << endl;
  o << decindent;
}

bool
Extern_RefWavefunction::sdref() const {
  // references based on Extern_RefWavefunction are single-determinant references IF their densities are idempotent
  return this->ordm_idempotent();
}

RefSymmSCMatrix
Extern_RefWavefunction::core_hamiltonian_for_basis(const Ref<GaussianBasisSet> &basis,
                                                   const Ref<GaussianBasisSet> &p_basis) {
  assert(p_basis == 0); // can only do nonrelativistic Hamiltonians now
  Ref<OrbitalSpaceRegistry> oreg = this->world()->tfactory()->orbital_registry();
  Ref<AOSpaceRegistry> aoreg = this->world()->tfactory()->ao_registry();
  const bool need_to_add_aospace_temporarily = !aoreg->key_exists(basis);
  if (need_to_add_aospace_temporarily)
    add_ao_space(basis, integral(), aoreg, oreg);

  const Ref<OrbitalSpace>& aox = aoreg->value(basis);
  const std::string nonrel_hkey =
    ParsedOneBodyIntKey::key(aox->id(), aox->id(),
                             std::string("H"));
  RefSCMatrix Hnr = world()->fockbuild_runtime()->get(nonrel_hkey);

  RefSymmSCMatrix Hnr_symm = Hnr.kit()->symmmatrix(Hnr.rowdim());
  Hnr_symm.assign(0.0);
  Hnr_symm.accumulate_symmetric_sum(Hnr);
  Hnr_symm.scale(0.5);

  if (need_to_add_aospace_temporarily)
    remove_ao_space(basis, aoreg, oreg);

  Ref<Integral> localints = integral()->clone();
  localints->set_basis(basis);
  Ref<PetiteList> pl = localints->petite_list();
  RefSymmSCMatrix Hnr_symm_so = pl->to_SO_basis(Hnr_symm);

  return Hnr_symm_so;
}

void
Extern_RefWavefunction::init(const RefSCMatrix& coefs,
                             const std::vector<unsigned int>& orbsyms,
                             std::vector<unsigned int> occpi,
                             std::vector<unsigned int> fzcpi,
                             std::vector<unsigned int> fzvpi,
                             std::vector<unsigned int> holepi,
                             std::vector<unsigned int> partpi)
{
  if (spinspaces_[Alpha].null()) {
    Extern_RefWavefunction* this_nonconst = const_cast<Extern_RefWavefunction*>(this);

    // make sure that FockBuildRuntime uses same density fitting info as this reference
    world()->fockbuild_runtime()->dfinfo(dfinfo());

    // make sure that FockBuildRuntime uses same densities as the reference wavefunction
    if(force_average_AB_rdm1_ == false) // the densites are in AO basis
        world()->fockbuild_runtime()->set_densities(ordm(Alpha), ordm(Beta));//here computes ordm
    else
    {
        RefSymmSCMatrix av_rdm = ordm(Alpha).copy();
        av_rdm.accumulate(ordm(Beta));
        av_rdm.scale(0.5);
        world()->fockbuild_runtime()->set_densities(av_rdm, av_rdm);
    }

    // make spaces
    init_spaces(coefs, orbsyms, occpi, fzcpi, fzvpi, holepi, partpi);
  }
}

void
Extern_RefWavefunction::init_spaces(const RefSCMatrix& coefs,
                                    const std::vector<unsigned int>& orbsyms,
                                    std::vector<unsigned int> occpi,
                                    std::vector<unsigned int> fzcpi,
                                    std::vector<unsigned int> fzvpi,
                                    std::vector<unsigned int> holepi,
                                    std::vector<unsigned int> partpi) {

  // compute orbspi
  const unsigned int nirreps = *std::max_element(orbsyms.begin(), orbsyms.end())
      + 1;
  std::vector<unsigned int> orbspi(nirreps);
  for (unsigned int irrep = 0; irrep < nirreps; ++irrep) {
    orbspi[irrep] = std::count(orbsyms.begin(), orbsyms.end(), irrep);
  }

  assert(occpi.size() == nirreps);
  assert(fzcpi.size() == nirreps);
  assert(fzvpi.size() == nirreps);

  // compute holepi and partpi, if needed
  if (holepi.empty()) {
    holepi.resize(nirreps);
    for (unsigned int irrep = 0; irrep < nirreps; ++irrep)
      holepi[irrep] = occpi[irrep] - fzcpi[irrep];
  }
  if (partpi.empty()) {
    partpi.resize(nirreps);
    for (unsigned int irrep = 0; irrep < nirreps; ++irrep)
      partpi[irrep] = orbspi[irrep] - occpi[irrep] - fzvpi[irrep];
  }

  assert(holepi.size() == nirreps);
  assert(partpi.size() == nirreps);

  const unsigned int nmo = coefs.coldim().n();
  const bool debug = false; // control debug printing
  const bool debug_version = true; // experiment with debug code
  RefSCMatrix C_ao;
  RefDiagSCMatrix evals;
  RefDiagSCMatrix corr_evals;
  // block orbitals by symmetry
  Ref<OrbitalSpace> pspace_ao;
  {
    const std::string unique_id = new_unique_key(
        this->world()->tfactory()->orbital_registry());
    const std::string id = ParsedOrbitalSpaceKey::key(unique_id, AnySpinCase1);
    const std::string name("MOs blocked by symmetry");
    evals = coefs.kit()->diagmatrix(coefs.coldim());
    RefDiagSCMatrix occnums = evals.clone();
    corr_evals = evals.clone(); // a helper mat to generate correlated orbs
    // will set eigenvalues as follows:
    // frozen occupied orbitals to -2.0
    // active occupied orbitals to -1.0
    // active virtual orbitals to 1.0
    // frozen virtual orbitals to 2.0
    // this will help to determine occupation numbers later
    occnums.assign(0.0);
    corr_evals.assign(2.0);
    evals.assign(2.0);
    if (debug) {
      evals->print("debug ref_wfn: print evals");
      occnums->print("debug ref_wfn: print occnums");
    }
    const unsigned int nirrep = occpi.size();
    int mocount = 0;
    for (int g = 0; g < nirrep; ++g) {
      for (int i = 0; i < fzcpi[g]; ++i) {
        evals.set_element(mocount + i, -2.0);
        corr_evals.set_element(mocount + i, -2.0);
        occnums.set_element(mocount + i, 1.0);
      }
      if (debug) {
        evals->print("debug ref_wfn: print evals");
        occnums->print("debug ref_wfn: print occnums");
      }
      for (int j = fzcpi[g]; j < occpi[g]; ++j) {
        evals.set_element(mocount + j, -1.0);
        occnums.set_element(mocount + j, 1.0);
        if (j - fzcpi[g] < holepi[g])
          corr_evals.set_element(mocount + j, -2.0);
      }
      if (debug) {
        evals->print("debug ref_wfn: print evals");
        occnums->print("debug ref_wfn: print occnums");
      }
      for (int k = occpi[g]; k < (orbspi[g] - fzvpi[g]); ++k) {
        evals.set_element(mocount + k, 1.0);
      }
      mocount += orbspi[g];
    }
    if (debug) {
      evals->print("debug ref_wfn: print evals");
      occnums->print("debug ref_wfn: print occnums");
    }
    pspace_ao = new OrderedOrbitalSpace<SymmetryMOOrder>(
        id, name, basis(), integral(), coefs, evals, occnums, orbsyms,
        SymmetryMOOrder(nirrep));
    if (debug) {
      pspace_ao->print_detail();
      pspace_ao->evals().print("debu ref_wfn: pspace_ao");
    }
  }

  if (debug_version) {
    C_ao = coefs;
  } else {
    C_ao = pspace_ao->coefs();
    evals = pspace_ao->evals();
  }
  if (debug) {
    evals.print("debug ref_wfn: pspace_ao");
  }

  // compute occupancies from evals (see the trick above)
  std::vector<double> occnums(nmo, 0.0);
  for (unsigned int i = 0; i < nmo; ++i) {
    occnums[i] = evals(i) > 0.0 ? 0.0 : 1.0;
  }

  // compute correlation-active orbitals
  std::vector<ParticleHoleOrbitalAttributes> actmask(nmo, ParticleHoleOrbitalAttributes::None);
  {
    size_t offset = 0;
    for(unsigned int irrep=0; irrep<nirreps; ++irrep) {
      std::fill(actmask.begin() + offset + fzcpi[irrep],
                actmask.begin() + offset + fzcpi[irrep] + holepi[irrep],
                ParticleHoleOrbitalAttributes::Hole);
      std::fill(actmask.begin() + offset + orbspi[irrep] - fzvpi[irrep] - partpi[irrep],
                actmask.begin() + offset + orbspi[irrep] - fzvpi[irrep],
                ParticleHoleOrbitalAttributes::Particle);

      offset+=orbspi[irrep];
    }
  }

  // transform orbitals to SO basis
  Ref<PetiteList> plist = integral()->petite_list();
  RefSCMatrix C_so = plist->evecs_to_SO_basis(C_ao);

  // compute overlap in SO basis
  RefSymmSCMatrix S_so = compute_onebody_matrix<&Integral::overlap>(plist);

  // transform to MO basis, verify orthonormality of MOs
  RefSymmSCMatrix S_mo = S_so.kit()->symmmatrix(C_so.coldim());
  S_mo.assign(0.0);
  S_mo->accumulate_transform(C_so, S_so, SCMatrix::TransposeTransform);
  //S_mo.print("MO metric");

  // compute "real" orbital eigenvalues as diagonal elements of the Fock matrix
#if 1
  {
    Ref<OrbitalSpaceRegistry> oreg =
        this->world()->tfactory()->orbital_registry();
    const bool need_to_add_pspace_temporarily = !oreg->key_exists(
        pspace_ao->id());
    if (need_to_add_pspace_temporarily)
      oreg->add(pspace_ao->id(), pspace_ao);
    Ref<AOSpaceRegistry> aoreg = this->world()->tfactory()->ao_registry();
    const bool need_to_add_aospace_temporarily = !aoreg->key_exists(basis());
    if (need_to_add_aospace_temporarily)
      add_ao_space(basis(), integral(), aoreg, oreg);

    const std::string fkey = ParsedOneBodyIntKey::key(pspace_ao->id(),
                                                      pspace_ao->id(),
                                                      std::string("F"));
    RefSCMatrix F = world()->fockbuild_runtime()->get(fkey);
    for (int o = 0; o < nmo; ++o)
      evals.set_element(o, F(o, o));
    //F.print("orbital Fock matrix");
    //evals.print("orbital eigenvalues");
    RefSymmSCMatrix P = world()->fockbuild_runtime()->P();
    RefSymmSCMatrix P_mo = P.kit()->symmmatrix(C_ao.coldim());
    P_mo.assign(0.0);
    P_mo.accumulate_transform(C_ao, P, SCMatrix::TransposeTransform);

    if (need_to_add_pspace_temporarily)
      oreg->remove(pspace_ao->id());
    if (need_to_add_aospace_temporarily)
      remove_ao_space(basis(), aoreg, oreg);

    //this->world()->obsolete();
  }
#endif

  Ref<OrbitalSpaceRegistry> oreg =
      this->world()->tfactory()->orbital_registry();
  spinspaces_[Alpha] = new PopulatedOrbitalSpace(oreg, AnySpinCase1, basis(),
                                                 integral(), C_ao, occnums,
                                                 actmask, evals);
  spinspaces_[Beta] = spinspaces_[Alpha];
#if 1
  spinspaces_[Alpha]->occ_sb()->print_detail();
  spinspaces_[Alpha]->occ_act_sb()->print_detail();
  spinspaces_[Alpha]->uocc_sb()->print_detail();
  spinspaces_[Alpha]->uocc_act_sb()->print_detail();
  spinspaces_[Alpha]->occ()->print_detail();
  spinspaces_[Alpha]->occ_act()->print_detail();
#endif

}

double
Extern_RefWavefunction::magnetic_moment() const {
  return magmom_;
}

Ref<DensityFittingInfo>
Extern_RefWavefunction::dfinfo() const {
  return use_world_dfinfo() ? const_cast<DensityFittingInfo*>(world()->tfactory()->df_info()) : 0;
  //return const_cast<DensityFittingInfo*>(world()->tfactory()->df_info());
}

///////////////////////////////////////////////////////////////////

#if !defined(HAVE_PSI3)
Ref<RefWavefunction>
RefWavefunctionFactory::make(const Ref<WavefunctionWorld> & world,
                                const Ref<Wavefunction> & ref,
                                bool spin_restricted,
                                unsigned int nfzc,
                                unsigned int nfzv,
                                Ref<OrbitalSpace> vir_space)
{
  { // OneBodyWavefunction
    Ref<OneBodyWavefunction> cast; cast << ref;
    if (cast.nonnull())
      return new SD_RefWavefunction(world, cast, spin_restricted, nfzc, nfzv, vir_space);
  }
  throw FeatureNotImplemented("this reference wavefunction cannot be used for R12 methods",
                              __FILE__, __LINE__);
}
#endif



