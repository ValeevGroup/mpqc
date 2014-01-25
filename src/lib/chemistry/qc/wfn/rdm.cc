//
// rdm.cc
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

#include <cassert>
#include <chemistry/qc/wfn/rdm.h>
#include <chemistry/qc/wfn/orbitalspace.h>
#include <math/mmisc/pairiter.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/basis/obintfactory.h>
#include <cassert>

using namespace sc;

namespace {
  /// compute MO density of a given spin
  RefSymmSCMatrix
  mo_density(const Ref<OneBodyWavefunction>& obwfn, SpinCase1 spin) {
    if ((spin == Beta || spin == AnySpinCase1) && !obwfn->spin_polarized())
      return mo_density(obwfn, Alpha);

    if (spin == AnySpinCase1 && obwfn->spin_polarized())
      ProgrammingError("asked for any spin density but the density is spin-polarized", __FILE__, __LINE__);

    RefSymmSCMatrix P_mo = obwfn->basis_matrixkit()->symmmatrix(obwfn->oso_dimension());
    P_mo.assign(0.0);
    const int nmo = P_mo.n();
    for(int mo=0; mo<nmo; ++mo)
      P_mo.set_element(mo, mo, (spin == Alpha) ? obwfn->alpha_occupation(mo)
                                               : obwfn->beta_occupation(mo) );

#if 0
      P_mo.print(prepend_spincase(spin,"OneBodyWavefunction MO density").c_str());
#endif

    return P_mo;
  }
}

//////////////////////

namespace sc {

template<>
size_t
RDM<One>::ndim(SpinCase1 spincase) const {
  return this->orbs(spincase)->rank();
}

template<>
Ref< RDMCumulant<One> >
RDM<One>::cumulant() const
{
  return new RDMCumulant<One>(const_cast<RDM<One>*>(this));
}

template<>
Ref< RDM<Zero> >
RDM<One>::rdm_m_1() const {
  throw ProgrammingError("RDM<One>::rdm_m_1() called",
                         __FILE__, __LINE__);
}

template<>
RefSymmSCMatrix
RDM<One>::scmat(SpinCase1 spin) const {
  if (scmat_[spin].nonnull()) return scmat_[spin];
  if (wfn_->spin_polarized() && spin == Beta)
    return scmat(Alpha);

  // need to transform density from AO basis to orbs basis
  // P' = C^t S P S C
  Ref<OrbitalSpace> orbs = this->orbs(spin);
  RefSCMatrix C = orbs->coefs();
  RefSymmSCMatrix P_ao = (spin == Alpha) ? wfn()->alpha_ao_density()
                                         : wfn()->beta_ao_density();

  Ref<PetiteList> plist = wfn()->integral()->petite_list();
  RefSymmSCMatrix S_so = compute_onebody_matrix<&Integral::overlap>(plist);
  RefSymmSCMatrix S_ao = plist->to_AO_basis(S_so);
  S_so = 0;
  RefSymmSCMatrix SPS_ao = P_ao.kit()->symmmatrix(S_ao.dim()); SPS_ao.assign(0.0);
  SPS_ao.accumulate_transform(S_ao, P_ao);

  scmat_[spin] = C.kit()->symmmatrix(C.coldim());
  scmat_[spin].assign(0.0);
  scmat_[spin].accumulate_transform(C, SPS_ao, SCMatrix::TransposeTransform);

  return scmat_[spin];
}

////////////////////////

template<>
RefSymmSCMatrix
RDMCumulant<One>::scmat(SpinCase1 spin) const {
  return density_->scmat(spin);
}

////////////////////////

template<>
size_t
RDM<Two>::ndim(SpinCase2 spincase) const {
  const SpinCase1 spin1 = case1(spincase);
  const SpinCase1 spin2 = case2(spincase);
  const int n1 = this->orbs(spin1)->rank();
  const int n2 = this->orbs(spin2)->rank();
  switch (spincase) {
    case AlphaAlpha:
    case BetaBeta:
      return n1 * (n1-1) / 2;
    case AlphaBeta:
      return n1 * n2;
    default:
      MPQC_ASSERT(false); // should not be reachable
  }
  MPQC_ASSERT(false);  // not reachable
  return 0;  // dummy return statement to pacify picky compilers
}

template<>
RefSymmSCMatrix
RDM<Two>::scmat(SpinCase2 spin) const {
  throw ProgrammingError("RDM<Two>::scmat is called",
                         __FILE__,
                         __LINE__);
}

template<>
Ref< RDM<One> >
RDM<Two>::rdm_m_1() const {
  throw ProgrammingError("RDM<Two>::rdm_m_1 is called",
                         __FILE__,
                         __LINE__);
}

template<>
Ref< RDMCumulant<Two> >
RDM<Two>::cumulant() const {
  throw ProgrammingError("RDM<Two>::cumulant is called",
                         __FILE__,
                         __LINE__);
}

////////////////////////

template<>
RefSymmSCMatrix
RDMCumulant<Two>::scmat(SpinCase2 spin) const {

  if (scmat_[spin].nonnull()) return scmat_[spin];
  if (!density_->wfn()->spin_polarized() && spin == BetaBeta)
    return scmat(AlphaAlpha);

  const SpinCase1 spin1 = case1(spin);
  const SpinCase1 spin2 = case2(spin);
  Ref< RDM<One> > opdm_obj = density_->rdm_m_1();
  const RefSymmSCMatrix opdm1 = opdm_obj->scmat(spin1);
  const RefSymmSCMatrix opdm2 = opdm_obj->scmat(spin2);
  const RefSymmSCMatrix tpdm = density_->scmat(spin);
  const int n = opdm1.n();
  RefSymmSCMatrix lambda = tpdm.copy();

  Ref<OrbitalSpace> orbs1 = opdm_obj->orbs(spin1);
  Ref<OrbitalSpace> orbs2 = opdm_obj->orbs(spin2);
  SpinMOPairIter pq_iter(orbs1->rank(), orbs2->rank(), spin);
  SpinMOPairIter rs_iter(orbs1->rank(), orbs2->rank(), spin);
  for(pq_iter.start(); int(pq_iter); pq_iter.next()) {
    const int p = pq_iter.i();
    const int q = pq_iter.j();
    const int pq = pq_iter.ij();
    for(rs_iter.start(); int(rs_iter); rs_iter.next()) {
      const int r = rs_iter.i();
      const int s = rs_iter.j();
      const int rs = rs_iter.ij();
      if (pq < rs) continue;
      lambda.accumulate_element(pq,rs,-opdm1.get_element(p,r)*opdm2.get_element(q,s));
      if(spin!=AlphaBeta) {
        lambda.accumulate_element(pq,rs,opdm1.get_element(q,r)*opdm2.get_element(p,s));
      }
    }
  }

  //lambda.print(prepend_spincase(spin,"2-RDM cumulant").c_str());

  scmat_[spin] = lambda;
  return scmat_[spin];
}

//////////////////////

template<>
size_t
SpinFreeRDM<One>::ndim() const {
  return this->orbs()->rank();
}

template<>
Ref< SpinFreeRDM<Zero> >
SpinFreeRDM<One>::rdm_m_1() const {
  throw ProgrammingError("SpinFreeRDM<One>::rdm_m_1() called",
                         __FILE__, __LINE__);
}

template<>
RefSymmSCMatrix
SpinFreeRDM<One>::scmat() const {
  if (scmat_.nonnull()) return scmat_;

  // need to transform density from AO basis to orbs basis
  // P' = C^t S P S C
  Ref<OrbitalSpace> orbs = this->orbs();
  RefSCMatrix C = orbs->coefs();
  RefSymmSCMatrix P_ao = wfn()->alpha_ao_density() + wfn()->beta_ao_density();

  Ref<PetiteList> plist = wfn()->integral()->petite_list();
  RefSymmSCMatrix S_so = compute_onebody_matrix<&Integral::overlap>(plist);
  RefSymmSCMatrix S_ao = plist->to_AO_basis(S_so);
  S_so = 0;
  RefSymmSCMatrix SPS_ao = P_ao.kit()->symmmatrix(S_ao.dim()); SPS_ao.assign(0.0);
  SPS_ao.accumulate_transform(S_ao, P_ao);

  scmat_ = C.kit()->symmmatrix(C.coldim());
  scmat_.assign(0.0);
  scmat_.accumulate_transform(C, SPS_ao, SCMatrix::TransposeTransform);

  return scmat_;
}

template<>
const Ref<DistArray4>&
SpinFreeRDM<One>::da4() const {
  throw ProgrammingError("SpinFreeRDM<One>::da4() is called",
                         __FILE__, __LINE__);
  return da4_;
}

////////////////////////

template<>
size_t
SpinFreeRDM<Two>::ndim() const {
  const int n = this->orbs()->rank();
  return n * n;
}

template<>
RefSymmSCMatrix
SpinFreeRDM<Two>::scmat() const {
  throw ProgrammingError("SpinFreeRDM<Two>::scmat is called",
                         __FILE__,
                         __LINE__);
  return scmat_; // unreachable
}

template<>
const Ref<DistArray4>&
SpinFreeRDM<Two>::da4() const {
  throw ProgrammingError("SpinFreeRDM<Two>::da4 is called",
                         __FILE__,
                         __LINE__);
  return da4_; // unreachable
}

template<>
Ref< SpinFreeRDM<One> >
SpinFreeRDM<Two>::rdm_m_1() const {
  throw ProgrammingError("SpinFreeRDM<Two>::da4 is called",
                         __FILE__,
                         __LINE__);
}

} // end of namespace sc

//////////////////////

ClassDesc
OBWfnRDMTwo::class_desc_(typeid(OBWfnRDMTwo),
                     "OBWfnRDMTwo",
                     1,               // version
                     "public RDM<Two>", // must match parent
                     0,               // change to create<OBWfnRDMTwo> if this class is DefaultConstructible
                     create<OBWfnRDMTwo>, // change to 0 if this class is not KeyValConstructible
                     create<OBWfnRDMTwo>  // change to 0 if this class is not StateInConstructible
                     );

OBWfnRDMTwo::OBWfnRDMTwo(const Ref<KeyVal>& kv) : RDM<Two>(kv) {
  wfn_ = require_dynamic_cast<OneBodyWavefunction*>(
        kv->describedclassvalue("wfn").pointer(),
        "OBWfnRDMTwo::OBWfnRDMTwo\n"
        );
}

OBWfnRDMTwo::OBWfnRDMTwo(StateIn& si) : RDM<Two>(si) {
  wfn_ << SavableState::restore_state(si);
  if (wfn_.null())
    throw ProgrammingError("failed constructor",__FILE__,__LINE__);
}

OBWfnRDMTwo::OBWfnRDMTwo(const Ref<OneBodyWavefunction>& wfn) : RDM<Two>(wfn.pointer()) {
}

OBWfnRDMTwo::~OBWfnRDMTwo() {
}

void
OBWfnRDMTwo::save_data_state(StateOut& so) {
  SavableState::save_state(wfn_.pointer(), so);
}

Ref<OBWfnRDMTwo::cumulant_type> sc::OBWfnRDMTwo::cumulant() const
{
  return new OBWfnRDMCumulantTwo(const_cast<OBWfnRDMTwo*>(this));
}

Ref<OBWfnRDMTwo::rdm_m_1_type> sc::OBWfnRDMTwo::rdm_m_1() const
{
  return new OBWfnRDMOne(this->wfn());
}

RefSymmSCMatrix sc::OBWfnRDMTwo::scmat(SpinCase2 spin) const {

  if (scmat_[spin].nonnull()) return scmat_[spin];
  if (!wfn_->spin_polarized() && spin == BetaBeta)
    return scmat(AlphaAlpha);

  Ref<SCMatrixKit> kit = SCMatrixKit::default_matrixkit();
  const SpinCase1 spin1 = case1(spin);
  const SpinCase1 spin2 = case2(spin);
  const RefSymmSCMatrix opdm1 = mo_density(wfn_, spin1);
  const RefSymmSCMatrix opdm2 = (wfn_->spin_polarized() && spin == AlphaBeta)
                                  ? mo_density(wfn_, spin2)
                                  : opdm1;
  const int n = opdm1.n();

  RefSCDimension dim = new SCDimension( spin == AlphaBeta ? n*n : n*(n-1)/2 );
  RefSymmSCMatrix tpdm = kit->symmmatrix(dim); tpdm.assign(0.0);

  int b12 = 0;
  for (int b1 = 0; b1 < n; ++b1) {
    const int b2fence = (spin == AlphaBeta) ? n : b1;
    for (int b2 = 0; b2 < b2fence; ++b2, ++b12) {

      int k12 = 0;
      for (int k1 = 0; k1 < n; ++k1) {
        const double gamma_b1_k1 = opdm1.get_element(b1, k1);
        const double gamma_b2_k1 = (spin != AlphaBeta) ? opdm1.get_element(b2,
                                                                           k1)
                                                       : 0.0;

        const int k2fence = (spin == AlphaBeta) ? n : k1;
        for (int k2 = 0; k2 < k2fence; ++k2, ++k12) {
          double value = gamma_b1_k1 * opdm2.get_element(b2, k2);
          if (spin != AlphaBeta)
            value -= opdm1.get_element(b1, k2) * gamma_b2_k1;
          tpdm.accumulate_element(b12, k12, value);
        }
      }
    }
  }
  scmat_[spin] = tpdm;
  return scmat_[spin];
}

namespace {
  Ref<OrbitalSpace> orbs_from_obwfn(const Ref<OneBodyWavefunction>& wfn,
                                    SpinCase1 spin) {
    const Ref<GaussianBasisSet> bs = wfn->basis();
    RefSCMatrix evecs_so = (spin == Alpha) ? wfn->alpha_eigenvectors()
                                           : wfn->beta_eigenvectors();
    const RefDiagSCMatrix evals = (spin == Alpha) ? wfn->alpha_eigenvalues()
                                                  : wfn->beta_eigenvalues();
    const Ref<Integral>& integral = wfn->integral();
    Ref<PetiteList> plist = integral->petite_list();
    const RefSCMatrix evecs_ao = plist->evecs_to_AO_basis(evecs_so);
    evecs_so = 0;

    const std::string prefix(to_string(spin));
    std::ostringstream oss;
    oss << prefix << " symmetry-blocked MOs";
    std::string id = ParsedOrbitalSpaceKey::key(std::string("p(sym)"), spin);
    Ref<OrbitalSpace> orbs = new OrbitalSpace(id, oss.str(), evecs_ao, bs, integral,
                                              evals, 0, 0, OrbitalSpace::symmetry);

    return orbs;
  }
}

Ref<OrbitalSpace>
sc::OBWfnRDMTwo::orbs(SpinCase1 spin) const {
  if (orbs_[spin].nonnull()) return orbs_[spin];
  if (wfn()->spin_polarized() && spin == Beta)
    return orbs(Alpha);

  orbs_[spin] = orbs_from_obwfn(wfn(), spin);
  return orbs_[spin];
}

/////////////////////

ClassDesc
OBWfnRDMCumulantTwo::class_desc_(typeid(OBWfnRDMCumulantTwo),
                   "OBWfnRDMCumulantTwo",
                   1,               // version
                   "public RDMCumulant<Two>", // must match parent
                   0,               // change to create<OBWfnRDMCumulantTwo> if this class is DefaultConstructible
                   0, // change to 0 if this class is not KeyValConstructible
                   create<OBWfnRDMCumulantTwo>  // change to 0 if this class is not StateInConstructible
    );

OBWfnRDMCumulantTwo::OBWfnRDMCumulantTwo(const Ref<OBWfnRDMTwo>& density) : density_(density), RDMCumulant<Two>(density) {
}

OBWfnRDMCumulantTwo::OBWfnRDMCumulantTwo(StateIn& si) : RDMCumulant<Two>(si) {
  density_ << SavableState::restore_state(si);
  if (density_.null())
    throw ProgrammingError("failed constructor",__FILE__,__LINE__);
}

OBWfnRDMCumulantTwo::~OBWfnRDMCumulantTwo() {
}

void
OBWfnRDMCumulantTwo::save_data_state(StateOut& so) {
  SavableState::save_state(density_.pointer(), so);
}

void sc::OBWfnRDMCumulantTwo::release_block(SpinCase2 spin, size_t bra, double *blk) const
{
  throw "not yet implemented";
}

const double *sc::OBWfnRDMCumulantTwo::obtain_block(SpinCase2 spin, size_t bra) const
{
  throw "not yet implemented";
}

void sc::OBWfnRDMCumulantTwo::compute()
{
  density_->compute();
}

RefSymmSCMatrix sc::OBWfnRDMCumulantTwo::scmat(SpinCase2 spin) const {

  if (scmat_[spin].nonnull()) return scmat_[spin];
  if (!density_->wfn()->spin_polarized() && spin == BetaBeta)
    return scmat(AlphaAlpha);

  const SpinCase1 spin1 = case1(spin);
  const SpinCase1 spin2 = case2(spin);
  const RefSymmSCMatrix opdm1 = mo_density(density_->wfn(), spin1);
  const RefSymmSCMatrix opdm2 = (density_->wfn()->spin_polarized() && spin == AlphaBeta)
                                  ? mo_density(density_->wfn(), spin2)
                                  : opdm1;
  const RefSymmSCMatrix tpdm = density_->scmat(spin);
  const int n = opdm1.n();
  RefSymmSCMatrix lambda = tpdm.copy();

  int b12 = 0;
  for(int b1=0; b1<n; ++b1) {
    const int b2fence = (spin == AlphaBeta) ? n : b1;
    for(int b2=0; b2< b2fence; ++b2, ++b12) {

      int k12 = 0;
      for(int k1=0; k1<n; ++k1) {
        const double gamma_b1_k1 = opdm1.get_element(b1,k1);
        const double gamma_b2_k1 = (spin != AlphaBeta) ? opdm1.get_element(b2,k1) : 0.0;

        const int k2fence = (spin == AlphaBeta) ? n : k1;
        for(int k2=0; k2< k2fence; ++k2, ++k12) {
          double value = - (gamma_b1_k1 * opdm2.get_element(b2,k2));
          if (spin != AlphaBeta) value += opdm1.get_element(b1,k2) * gamma_b2_k1;
          lambda.accumulate_element(b12,k12,value);
        }
      }
    }
  }

  scmat_[spin] = lambda;
  return scmat_[spin];
}

#if 0
/////////////////////

ClassDesc
WfnRDMOne::class_desc_(typeid(WfnRDMOne),
                       "WfnRDMOne",
                       1,               // version
                       "public RDM<One>", // must match parent
                       0, 0, 0
                     );

WfnRDMOne::WfnRDMOne(const Ref<KeyVal>& kv) : RDM<One>(kv) {
  wfn_ = require_dynamic_cast<Wavefunction*>(
        kv->describedclassvalue("wfn").pointer(),
        "WfnRDMOne::WfnRDMOne\n"
        );
}

WfnRDMOne::WfnRDMOne(StateIn& si) : RDM<One>(si) {
  wfn_ << SavableState::restore_state(si);
  if (wfn_.null())
    throw ProgrammingError("failed constructor",__FILE__,__LINE__);
}

WfnRDMOne::WfnRDMOne(const Ref<Wavefunction>& wfn) : RDM<One>(wfn) {
}

WfnRDMOne::~WfnRDMOne() {
}

void
WfnRDMOne::save_data_state(StateOut& so) {
  SavableState::save_state(wfn_.pointer(), so);
}

Ref<WfnRDMOne::cumulant_type> sc::WfnRDMOne::cumulant() const
{
  return new WfnRDMCumulantOne(const_cast<WfnRDMOne*>(this));
}

Ref< RDM<Zero> > sc::WfnRDMOne::rdm_m_1() const {
  throw ProgrammingError("RDM<One>::rdm_m_1() called",
                         __FILE__, __LINE__);
}

void sc::WfnRDMOne::release_block(SpinCase1 spin, size_t bra, double *blk) const
{
  throw "not yet implemented";
}

const double *sc::WfnRDMOne::obtain_block(SpinCase1 spin, size_t bra) const
{
  throw "not yet implemented";
}

void sc::WfnRDMOne::compute()
{
  // force computation of the value
  const double value = wfn_->value();
}

size_t sc::WfnRDMOne::ndim(SpinCase1 spincase) const
{
  return this->orbs(spincase)->rank();
}

/////////////////////

ClassDesc
WfnRDMCumulantOne::class_desc_(typeid(WfnRDMCumulantOne),
                   "WfnRDMCumulantOne",
                   1,               // version
                   "public RDMCumulant<One>", // must match parent
                   0,               // change to create<WfnRDMCumulantOne> if this class is DefaultConstructible
                   0, // change to 0 if this class is not KeyValConstructible
                   create<WfnRDMCumulantOne>  // change to 0 if this class is not StateInConstructible
    );

WfnRDMCumulantOne::WfnRDMCumulantOne(const Ref<WfnRDMOne>& density) :
    density_(density), RDMCumulant<One>(density) {
}

WfnRDMCumulantOne::WfnRDMCumulantOne(StateIn& si) : RDMCumulant<One>(si) {
  density_ << SavableState::restore_state(si);
  if (density_.null())
    throw ProgrammingError("failed constructor",__FILE__,__LINE__);
}

WfnRDMCumulantOne::~WfnRDMCumulantOne() {
}

void
WfnRDMCumulantOne::save_data_state(StateOut& so) {
  SavableState::save_state(density_.pointer(), so);
}

void sc::WfnRDMCumulantOne::release_block(SpinCase1 spin, size_t bra, double *blk) const
{
  throw "not yet implemented";
}

const double *sc::WfnRDMCumulantOne::obtain_block(SpinCase1 spin, size_t bra) const
{
  throw "not yet implemented";
}

void sc::WfnRDMCumulantOne::compute()
{
  density_->compute();
}

RefSymmSCMatrix sc::WfnRDMCumulantOne::scmat(SpinCase1 spin) const {
  return density_->scmat(spin);
}
#endif
/////////////////////

ClassDesc
OBWfnRDMOne::class_desc_(typeid(OBWfnRDMOne),
                       "OBWfnRDMOne",
                       1,               // version
                       "public RDM<One>", // must match parent
                       0, create<OBWfnRDMOne>, create<OBWfnRDMOne>
                     );

OBWfnRDMOne::OBWfnRDMOne(const Ref<KeyVal>& kv) : RDM<One>(kv) {
  wfn_ = require_dynamic_cast<OneBodyWavefunction*>(
        kv->describedclassvalue("wfn").pointer(),
        "OBWfnRDMOne::OBWfnRDMOne\n"
        );
}

OBWfnRDMOne::OBWfnRDMOne(StateIn& si) : RDM<One>(si) {
  wfn_ << SavableState::restore_state(si);
  if (wfn_.null())
    throw ProgrammingError("failed constructor",__FILE__,__LINE__);
}

OBWfnRDMOne::OBWfnRDMOne(const Ref<OneBodyWavefunction>& wfn) : RDM<One>(wfn.pointer()), wfn_(wfn) {
}

OBWfnRDMOne::~OBWfnRDMOne() {
}

void
OBWfnRDMOne::save_data_state(StateOut& so) {
  SavableState::save_state(wfn_.pointer(), so);
}

Ref<OrbitalSpace>
sc::OBWfnRDMOne::orbs(SpinCase1 spin) const {
  if (orbs_[spin].nonnull()) return orbs_[spin];
  if (wfn()->spin_polarized() && spin == Beta)
    return orbs(Alpha);

  orbs_[spin] = orbs_from_obwfn(wfn(), spin);
  return orbs_[spin];
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
