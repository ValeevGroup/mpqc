//
// r12int_eval.cc
//
// Copyright (C) 2004 Edward Valeev
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
#include <util/misc/formio.h>
#include <util/misc/scexception.h>
#include <util/ref/ref.h>
#include <util/state/state_bin.h>
#include <math/scmat/local.h>
#include <chemistry/qc/wfn/wfn.h>
#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/mbptr12/r12wfnworld.h>
#include <math/mmisc/pairiter.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/lcao/transform_factory.h>
#include <chemistry/qc/lcao/utils.h>
#include <chemistry/qc/mbptr12/r12_amps.h>
#include <chemistry/qc/mbptr12/compute_tbint_tensor.h>
#include <chemistry/qc/mbptr12/contract_tbint_tensor.h>
#include <chemistry/qc/mbptr12/container.h>
#include <chemistry/qc/mbptr12/creator.h>
#include <chemistry/qc/mbptr12/debug.h>
#include <chemistry/qc/wfn/orbitalspace_utils.h>
#include <chemistry/qc/basis/obintfactory.h>

using namespace std;
using namespace sc;

#define INCLUDE_EBC_CODE 1
#define INCLUDE_COUPLING_CODE 1
#define INCLUDE_GBC_CODE 1

inline int max(int a,int b) { return (a > b) ? a : b;}

/*-----------------
  R12IntEval
 -----------------*/
static ClassDesc R12IntEval_cd(
  typeid(R12IntEval),"R12IntEval",5,"virtual public SavableState",
  0, 0, 0);

R12IntEval::R12IntEval(const Ref<R12WavefunctionWorld>& r12w) :
  r12world_(r12w), evaluated_(false), debug_(0), emp2_obs_singles_(1.0), emp2_cabs_singles_(1.0)
{
  int naocc_a, naocc_b;
  int navir_a, navir_b;
  int nall_a, nall_b;
  Ref<RefWavefunction> refinfo = r12world()->refwfn();

  const bool nalpha_eq_nbeta = r12world()->spinadapted() || (spin_polarized() == false && r12world()->spinadapted() == false);

  if (nalpha_eq_nbeta) {
    const int nocc_act = refinfo->occ_act()->rank(); // this will lead to RASCI computation and orbital space initilization
    const int nvir_act = refinfo->uocc_act()->rank();
    const int nall = refinfo->orbs()->rank();
    naocc_a = naocc_b = nocc_act;
    navir_a = navir_b = nvir_act;
    nall_a = nall_b = nall;
  }
  else {
    naocc_a = occ_act(Alpha)->rank();
    naocc_b = occ_act(Beta)->rank();
    navir_a = vir_act(Alpha)->rank();
    navir_b = vir_act(Beta)->rank();
    nall_a = r12world()->refwfn()->orbs(Alpha)->rank();
    nall_b = r12world()->refwfn()->orbs(Beta)->rank();
  }

  dim_oo_[AlphaAlpha] = new SCDimension((naocc_a*(naocc_a-1))/2);
  dim_vv_[AlphaAlpha] = new SCDimension((navir_a*(navir_a-1))/2);
  dim_aa_[AlphaAlpha] = new SCDimension((nall_a*(nall_a-1))/2);
  dim_oo_[AlphaBeta] = new SCDimension(naocc_a*naocc_b);
  dim_vv_[AlphaBeta] = new SCDimension(navir_a*navir_b);
  dim_aa_[AlphaBeta] = new SCDimension(nall_a*nall_b);
  dim_oo_[BetaBeta] = new SCDimension((naocc_b*(naocc_b-1))/2);
  dim_vv_[BetaBeta] = new SCDimension((navir_b*(navir_b-1))/2);
  dim_aa_[BetaBeta] = new SCDimension((nall_b*(nall_b-1))/2);

  Ref<R12Technology> r12tech = r12world()->r12tech();
  switch(r12tech->ansatz()->orbital_product_gg()) {
    case R12Technology::OrbProdgg_ij:
      dim_gg_[AlphaAlpha] = dim_oo_[AlphaAlpha];
      dim_gg_[AlphaBeta] = dim_oo_[AlphaBeta];
      dim_gg_[BetaBeta] = dim_oo_[BetaBeta];
      break;
    case R12Technology::OrbProdgg_pq:
      {
        const unsigned int norbs_a = this->orbs(Alpha)->rank();
        const unsigned int norbs_b = this->orbs(Beta)->rank();
        dim_gg_[AlphaAlpha] = new SCDimension((norbs_a*(norbs_a-1))/2);
        dim_gg_[AlphaBeta] = new SCDimension(norbs_a*norbs_b);
        dim_gg_[BetaBeta] = new SCDimension((norbs_b*(norbs_b-1))/2);
      }
      break;
    default:
        throw ProgrammingError("R12IntEval::R12IntEval -- invalid orbital_product_gg for the R12 ansatz",__FILE__,__LINE__);
  }
  switch(r12tech->ansatz()->orbital_product_GG()) {
  case R12Technology::OrbProdGG_ij:
      dim_GG_[AlphaAlpha] = dim_oo_[AlphaAlpha];
      dim_GG_[AlphaBeta] = dim_oo_[AlphaBeta];
      dim_GG_[BetaBeta] = dim_oo_[BetaBeta];
      break;
  case R12Technology::OrbProdGG_pq:
      {
        const unsigned int norbs_a = this->orbs(Alpha)->rank();
        const unsigned int norbs_b = this->orbs(Beta)->rank();
        dim_GG_[AlphaAlpha] = new SCDimension((norbs_a*(norbs_a-1))/2);
        dim_GG_[AlphaBeta] = new SCDimension(norbs_a*norbs_b);
        dim_GG_[BetaBeta] = new SCDimension((norbs_b*(norbs_b-1))/2);
      }
      break;
  default:
      throw ProgrammingError("R12IntEval::R12IntEval -- invalid orbital_product_GG for the R12 ansatz",__FILE__,__LINE__);
  }
  for(int s=0; s<NSpinCases2; s++) {
    dim_f12_[s] = new SCDimension(corrfactor()->nfunctions()*dim_GG_[s].n());
  }

  Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();
  for(int s=0; s<NSpinCases2; s++) {
    // if spin-free, calculation will only need to evaluate one type; but since many methods need to access all three (at least 2 spincases), we
    // still create AA and AB case
    // if spin-orbital, and spin-polarized allocate all three
    if ((r12world()->spinadapted() == true && s != BetaBeta) ||
        (r12world()->spinadapted() == false && spin_polarized() == false && s != BetaBeta) ||
        (r12world()->spinadapted() == false && spin_polarized() == true)
       ) {
      V_[s] = local_matrix_kit->matrix(dim_f12_[s],dim_gg_[s]);
      X_[s] = local_matrix_kit->matrix(dim_f12_[s],dim_f12_[s]);
      B_[s] = local_matrix_kit->matrix(dim_f12_[s],dim_f12_[s]);
      if (stdapprox() == R12Technology::StdApprox_B)
        BB_[s] = local_matrix_kit->matrix(dim_f12_[s],dim_f12_[s]);
      emp2pair_[s] = local_matrix_kit->vector(dim_gg_[s]);
      const SpinCase2 pairspin = static_cast<SpinCase2>(s);
      cuspconsistentgeminalcoefficient_[s] =
          new CuspConsistentGeminalCoefficient(pairspin);
    }
  }
  // if spin-orbital, but non-spin-polarized allocate AlphaBeta and AlphaAlpha (AlphaAlpha equiv to BetaBeta)
  if (spin_polarized() == false or r12world()->spinadapted()) {
    V_[BetaBeta] = V_[AlphaAlpha];
    X_[BetaBeta] = X_[AlphaAlpha];
    B_[BetaBeta] = B_[AlphaAlpha];
    BB_[BetaBeta] = BB_[AlphaAlpha];
    emp2pair_[BetaBeta] = emp2pair_[AlphaAlpha];
    cuspconsistentgeminalcoefficient_[BetaBeta] = cuspconsistentgeminalcoefficient_[AlphaAlpha];
  }

#if 0
  // WARNING Can use DistArray4_MemoryGrp only for MP2 computations
  bool mp2_only = false;
  {
    Ref<R12Technology::NullCorrelationFactor> nullptr; nullptr << corrfactor();
    if (nullptr)
      mp2_only = true;
  }
  if (r12world()->ints_method() != R12WavefunctionWorld::StoreMethod::posix && !mp2_only)
    throw InputError("R12IntEval::R12IntEval() -- the only supported storage method is posix");
#endif

  ordm_[Alpha] = ordm_[Beta] = 0;
}

R12IntEval::R12IntEval(StateIn& si) : SavableState(si)
{
  r12world_ << SavableState::restore_state(si);

  Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();
  for(int s=0; s<NSpinCases2; s++) {
    dim_oo_[s] << SavableState::restore_state(si);
    dim_vv_[s] << SavableState::restore_state(si);
    dim_f12_[s] << SavableState::restore_state(si);
    if (si.version(::class_desc<R12IntEval>()) >= 2) {
      dim_GG_[s] << SavableState::restore_state(si);
      dim_gg_[s] << SavableState::restore_state(si);
    }
    bool have_A;
    if (si.version(::class_desc<R12IntEval>()) >= 5) {
      si.get(have_A);
    }
    else {
      have_A = coupling() || ebc();
    }

    if (!(spin_polarized() && s == BetaBeta)) {
      V_[s] = local_matrix_kit->matrix(dim_f12_[s],dim_gg_[s]);
      X_[s] = local_matrix_kit->matrix(dim_f12_[s],dim_f12_[s]);
      B_[s] = local_matrix_kit->matrix(dim_f12_[s],dim_f12_[s]);
      if (stdapprox() == R12Technology::StdApprox_B)
        BB_[s] = local_matrix_kit->matrix(dim_f12_[s],dim_f12_[s]);
      if (have_A) {
        A_[s] = local_matrix_kit->matrix(dim_f12_[s],dim_vv_[s]);
      }
      emp2pair_[s] = local_matrix_kit->vector(dim_vv_[s]);
      const SpinCase2 pairspin = static_cast<SpinCase2>(s);
      cuspconsistentgeminalcoefficient_[s] = new CuspConsistentGeminalCoefficient(pairspin);

      V_[s].restore(si);
      X_[s].restore(si);
      B_[s].restore(si);
      BB_[s].restore(si);
      A_[s].restore(si);
      emp2pair_[s].restore(si);
    }
    else {
      V_[BetaBeta] = V_[AlphaAlpha];
      X_[BetaBeta] = X_[AlphaAlpha];
      B_[BetaBeta] = B_[AlphaAlpha];
      BB_[BetaBeta] = BB_[AlphaAlpha];
      A_[BetaBeta] = A_[AlphaAlpha];
      emp2pair_[BetaBeta] = emp2pair_[AlphaAlpha];
      cuspconsistentgeminalcoefficient_[BetaBeta] = cuspconsistentgeminalcoefficient_[AlphaAlpha];
    }
  }

  int evaluated; si.get(evaluated); evaluated_ = (bool) evaluated;
  si.get(debug_);
  if (si.version(::class_desc<R12IntEval>()) >= 3)
    si.get(emp2_obs_singles_);
  if (si.version(::class_desc<R12IntEval>()) >= 4)
    si.get(emp2_cabs_singles_);
}

R12IntEval::~R12IntEval()
{
}

void
R12IntEval::save_data_state(StateOut& so)
{
  SavableState::save_state(r12world_.pointer(),so);

  for(int s=0; s<NSpinCases2; s++) {
    SavableState::save_state(dim_oo_[s].pointer(),so);
    SavableState::save_state(dim_vv_[s].pointer(),so);
    SavableState::save_state(dim_f12_[s].pointer(),so);
    SavableState::save_state(dim_GG_[s].pointer(),so);
    SavableState::save_state(dim_gg_[s].pointer(),so);
    const bool have_A = A_[AlphaBeta];
    so.put(have_A);
    if (!(spin_polarized() && s == BetaBeta)) {
      V_[s].save(so);
      X_[s].save(so);
      B_[s].save(so);
      BB_[s].save(so);
      if (have_A) A_[s].save(so);
      emp2pair_[s].save(so);
    }
  }

  so.put((int)evaluated_);
  so.put(debug_);
  so.put(emp2_obs_singles_);
  so.put(emp2_cabs_singles_);
}

void
R12IntEval::obsolete()
{
  evaluated_ = false;

  emp2_obs_singles_ = 1.0;
  emp2_cabs_singles_ = 1.0;

  init_intermeds_();
}

const RefSCMatrix&
R12IntEval::V(SpinCase2 S) {
  compute();
  if (!spin_polarized() && (S == AlphaAlpha || S == BetaBeta))
    antisymmetrize(V_[AlphaAlpha],V_[AlphaBeta],
                   GGspace(Alpha),
                   ggspace(Alpha));
  return V_[S];
}


const RefSCMatrix&
R12IntEval::V() {
  MPQC_ASSERT(r12world()->spinadapted() == true);
  compute();
  return V_[AlphaBeta];
}


double R12IntEval::C_CuspConsistent(int i,int j,int k,int l,SpinCase2 pairspin) {
  return(cuspconsistentgeminalcoefficient_[pairspin]->C(i,j,k,l));
}

namespace {
  /// Returns the lower triangle of the matrix B (which should be symmetric)
  RefSymmSCMatrix to_lower_triangle(const RefSCMatrix& B) {
    RefSymmSCMatrix Bs = B.kit()->symmmatrix(B.rowdim());
    int n = B.nrow();
    double* b = new double[n*n];
    B.convert(b);
    const double* b_ptr = b;
    for(int i=0; i<n; i++, b_ptr += i)
      for(int j=i; j<n; j++, b_ptr++)
        Bs.set_element(i,j,*b_ptr);
    delete[] b;
    return Bs;
  }
}

RefSymmSCMatrix
R12IntEval::X(SpinCase2 S) {
  compute();
  if (!spin_polarized() && (S == AlphaAlpha || S == BetaBeta))
    antisymmetrize(X_[AlphaAlpha],X_[AlphaBeta],
                   GGspace(Alpha),
                   GGspace(Alpha));
  return to_lower_triangle(X_[S]);
}

RefSymmSCMatrix
R12IntEval::X() {
  MPQC_ASSERT(r12world()->spinadapted() == true);
  compute();
  return to_lower_triangle(X_[AlphaBeta]);
}

RefSymmSCMatrix
R12IntEval::B(SpinCase2 S) {
  compute();
  if (!spin_polarized() && (S == AlphaAlpha || S == BetaBeta))
    antisymmetrize(B_[AlphaAlpha],B_[AlphaBeta],
                   GGspace(Alpha),
                   GGspace(Alpha));
  return to_lower_triangle(B_[S]);
}

RefSymmSCMatrix
R12IntEval::B() {
  MPQC_ASSERT(r12world()->spinadapted() == true);
  compute();
  return to_lower_triangle(B_[AlphaBeta]);
}
RefSymmSCMatrix
R12IntEval::BB(SpinCase2 S) {
  if (stdapprox() != R12Technology::StdApprox_B)
    throw ProgrammingError("R12IntEval::BB() -- called but standard approximation is not B",__FILE__,__LINE__);
  compute();
  if (!spin_polarized() && (S == AlphaAlpha || S == BetaBeta))
    antisymmetrize(BB_[AlphaAlpha],BB_[AlphaBeta],
                   GGspace(Alpha),
                   GGspace(Alpha));
  return to_lower_triangle(BB_[S]);
}

const RefSCMatrix&
R12IntEval::A(SpinCase2 S) {
  if (A_[AlphaBeta].null())
    compute_A();
  if (!spin_polarized() && (S == AlphaAlpha || S == BetaBeta)) {
    if (A_[AlphaAlpha].null()) {
      Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();
      A_[AlphaAlpha] = local_matrix_kit->matrix(dim_f12_[AlphaAlpha],dim_vv_[AlphaAlpha]);
      antisymmetrize(A_[AlphaAlpha],A_[AlphaBeta],
                     GGspace(Alpha),
                     vir_act(Alpha));
    }
  }
  return A_[S];
}

const RefSCMatrix&
R12IntEval::T2(SpinCase2 S) {
  compute();
  return amps()->T2(S);
}

const RefSCMatrix&
R12IntEval::F12(SpinCase2 S) {
  compute();
  return amps()->Fvv(S);
}

Ref<R12Amplitudes>
R12IntEval::amps()
{
  if (Amps_.null())  // using this assumes this is not called from a constructor
    Amps_ = new R12Amplitudes(this);
  return Amps_;
}

const RefSCVector&
R12IntEval::compute_emp2(SpinCase2 spincase2)
{
  const int s = static_cast<int>(spincase2);
  // for single-determinant references:
  // 1) compute MP2 energies
  // 2) optional singles MP2 contribution for non-Brillouin references
  if (r12world()->sdref()) {
    const bool obs_eq_vbs = r12world()->obs_eq_vbs();
    const bool obs_eq_ribs = r12world()->obs_eq_ribs();
    Ref<R12Technology::NullCorrelationFactor> nocorrptr; nocorrptr << corrfactor();

    const SpinCase1 spin1 = case1(spincase2);
    const SpinCase1 spin2 = case2(spincase2);
    if (dim_oo(spincase2).n() == 0)
      return emp2pair_[s];
    Ref<OrbitalSpace> occ1_act = occ_act(spin1);
    Ref<OrbitalSpace> occ2_act = occ_act(spin2);
    Ref<OrbitalSpace> vir1_act = vir_act(spin1);
    Ref<OrbitalSpace> vir2_act = vir_act(spin2);

    std::string tform_key;
    // If VBS==OBS and this is not a pure MP2 calculation then this tform should be available
    if (obs_eq_vbs && nocorrptr.null()) {
      R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(), occ1_act,
                                               this->orbs(spin1), occ2_act,
                                               this->orbs(spin2), corrfactor(),
                                               true);
      tform_key = tformkey_creator();
    } else if (!obs_eq_vbs && nocorrptr.null()) { // MP2-R12 and VBS != OBS -- this transform will be available
      R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(), occ1_act,
                                               vir1_act, occ2_act, vir2_act,
                                               corrfactor(), true);
      tform_key = tformkey_creator();
    } else { // pure MP2 -- manually construct transform
      const std::string
          descr_key =
              moints_runtime4()->descr_key(
                                           new TwoBodyIntDescrERI(
                                                                  r12world()->integral()));
      const std::string layout_key = std::string(TwoBodyIntLayout::b1b2_k1k2);
      tform_key = ParsedTwoBodyFourCenterIntKey::key(occ1_act->id(),
                                                     occ2_act->id(),
                                                     vir1_act->id(),
                                                     vir2_act->id(), descr_key,
                                                     layout_key);
    }

    compute_mp2_pair_energies_(emp2pair_[s], spincase2, occ1_act, vir1_act,
                               occ2_act, vir2_act, tform_key);

    globally_sum_scvector_(emp2pair_[s], true);

  } // single-determinant references -> compute MP2 energies
  return emp2pair_[s];
}

double
R12IntEval::emp2_obs_singles()
{
  // compute OBS singles contribution to the MP2 energy if non-Brillouin reference is used
  if (!this->bc()) {
    const bool obs_singles = true;
    if (emp2_obs_singles_ == 1.0)
      emp2_obs_singles_ = compute_emp2_obs_singles(obs_singles);
  }
  return emp2_obs_singles_;
}

double
R12IntEval::emp2_cabs_singles(bool vir_cabs_coupling)
{
  if (!r12world()->sdref())
    ExEnv::out0() << indent
                  << "WARNING: CABS singles correction is not implemented for multiconfiguration references"
                  << std::endl;
  if (emp2_cabs_singles_ == 1.0)
    emp2_cabs_singles_ = compute_emp2_cabs_singles_noncanonical(vir_cabs_coupling);
  return emp2_cabs_singles_;
}

double
R12IntEval::emp2_cabs_singles(const RefSCMatrix& T1_ia_alpha,
                              const RefSCMatrix& T1_ia_beta) {
  MPQC_ASSERT(r12world()->sdref() == true);
  if (emp2_cabs_singles_ == 1.0)
    emp2_cabs_singles_ = compute_emp2_cabs_singles_noncanonical_ccsd(T1_ia_alpha,T1_ia_beta);
  return emp2_cabs_singles_;
}

const RefSCMatrix&
R12IntEval::T1_cabs(SpinCase1 spin) const {
  MPQC_ASSERT(r12world()->sdref() == true);
  MPQC_ASSERT(emp2_cabs_singles_ != 1.0);
  return T1_cabs_[spin];
}

const RefSCVector&
R12IntEval::emp2(SpinCase2 S)
{
  if (!spin_polarized() && S == BetaBeta)
    return compute_emp2(AlphaAlpha);
  else
    return compute_emp2(S);
}

void
R12IntEval::checkpoint_() const
{
  int me = r12world()->world()->msg()->me();
  Wavefunction* wfn = r12world()->wfn();

  if (me == 0 && wfn->if_to_checkpoint()) {
    StateOutBin stateout(wfn->checkpoint_file());
    SavableState::save_state(wfn,stateout);
    ExEnv::out0() << indent << "Checkpointed Wavefunction" << endl;
  }
}

void
R12IntEval::init_intermeds_()
{
  for(int s=0; s<NSpinCases2; s++) {
    V_[s].assign(0.0);
    X_[s].assign(0.0);
    B_[s].assign(0.0);
    if (BB_[s]) BB_[s].assign(0.0);
    emp2pair_[s].assign(0.0);
    if (A_[s]) A_[s].assign(0.0);
  }

  // nothing to do if no explicit correlation
  Ref<R12Technology::NullCorrelationFactor> no12ptr; no12ptr << corrfactor();
  if (no12ptr)
    return;

  Ref<R12Technology::G12CorrelationFactor> g12ptr; g12ptr << corrfactor();
  Ref<R12Technology::G12NCCorrelationFactor> g12ncptr; g12ncptr << corrfactor();
  Ref<R12Technology::R12CorrelationFactor> r12ptr; r12ptr << corrfactor();
  if (r12ptr) {
    init_intermeds_r12_();
  }
  else if (g12ptr || g12ncptr) {
    init_intermeds_g12_();
  }

  // add in the double-commutator relativistic contribution, if needed
  // can only use 1 Gaussian-based correlation factor (but pq ansatz should work!)
  if (this->dk() > 0)
    compute_B_DKH_();

}

void
R12IntEval::init_intermeds_r12_()
{
  for(int s=0; s<nspincases2(); s++) {

    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
    SpinCase1 spin1 = case1(spincase2);
    SpinCase1 spin2 = case2(spincase2);
	const Ref<OrbitalSpace>& GGspace1 = GGspace(case1(spincase2));
	const Ref<OrbitalSpace>& GGspace2 = GGspace(case2(spincase2));
    const Ref<OrbitalSpace>& gg1space = ggspace(spin1);
    const Ref<OrbitalSpace>& gg2space = ggspace(spin2);
    const Ref<OrbitalSpace>& GG1space = GGspace(spin1);
    const Ref<OrbitalSpace>& GG2space = GGspace(spin2);

	// compute identity operator in xc.pair/act.occ.pair basis
	RefSCMatrix I = compute_I_(GG1space,GG2space,gg1space,gg2space);
	if (spincase2 == AlphaBeta) {
	  V_[s].accumulate(I);
	}
	else {
	  antisymmetrize(V_[s],I,GG1space,gg1space);
	}

	if (r12world()->world()->msg()->me() == 0)
	  B_[s]->unit();
  }
  r2_contrib_to_X_new_();
}

/// Compute <space1 space2|space3 space4>
RefSCMatrix
R12IntEval::compute_I_(const Ref<OrbitalSpace>& space1,
		       const Ref<OrbitalSpace>& space2,
		       const Ref<OrbitalSpace>& space3,
		       const Ref<OrbitalSpace>& space4)
{
  /*--------------------------
    Compute overlap integrals
   --------------------------*/
  RefSCMatrix S_13 = compute_overlap_ints(space1, space3);
  RefSCMatrix S_24;
  if (space1 == space2 && space3 == space4) {
    S_24 = S_13;
  }
  else {
    S_24 = compute_overlap_ints(space2, space4);
  }
  const int nproc = r12world()->world()->msg()->n();
  const int me = r12world()->world()->msg()->me();

  const int n1 = space1->rank();
  const int n2 = space2->rank();
  const int n3 = space3->rank();
  const int n4 = space4->rank();
  const int n12 = n1*n2;
  const int n34 = n3*n4;
  const int n1234 = n12*n34;
  double* I_array = new double[n1234];
  memset(I_array,0,n1234*sizeof(double));

  int ij = 0;
  double* ijkl_ptr = I_array;
  for(int i=0; i<n1; i++)
    for(int j=0; j<n2; j++, ij++) {

    int ij_proc = ij%nproc;
    if (ij_proc != me) {
      ijkl_ptr += n34;
      continue;
    }

    int kl=0;
    for(int k=0; k<n3; k++)
      for(int l=0; l<n4; l++, kl++, ijkl_ptr++) {

        double S_ik = S_13.get_element(i,k);
        double S_jl = S_24.get_element(j,l);

        double I_ijkl = S_ik * S_jl;
        *ijkl_ptr = I_ijkl;
      }
    }

  r12world()->world()->msg()->sum(I_array,n1234);

  RefSCDimension dim_ij = new SCDimension(n12);
  RefSCDimension dim_kl = new SCDimension(n34);

  Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();
  RefSCMatrix I = local_matrix_kit->matrix(dim_ij, dim_kl);
  I.assign(I_array);
  delete[] I_array;

  // Only task 0 needs I
  if (me != 0)
    I.assign(0.0);

  return I;
}

/// Compute <space1 space2|r_{12}^2|space3 space4>
RefSCMatrix
R12IntEval::compute_r2_(const Ref<OrbitalSpace>& space1,
                        const Ref<OrbitalSpace>& space2,
                        const Ref<OrbitalSpace>& space3,
                        const Ref<OrbitalSpace>& space4)
{
  /*-----------------------------------------------------
    Compute overlap, dipole, quadrupole moment integrals
   -----------------------------------------------------*/
  RefSCMatrix S_13, MX_13, MY_13, MZ_13, MXX_13, MYY_13, MZZ_13, MXY_13, MXZ_13, MYZ_13;
  compute_multipole_ints(space1, space3,
                         MX_13, MY_13, MZ_13,
                         MXX_13, MYY_13, MZZ_13,
                         MXY_13, MXZ_13, MYZ_13);
  S_13 = compute_overlap_ints(space1, space3);

  RefSCMatrix S_24, MX_24, MY_24, MZ_24, MXX_24, MYY_24, MZZ_24, MXY_24, MXZ_24, MYZ_24;
  if (space1 == space2 && space3 == space4) {
    S_24 = S_13;
    MX_24 = MX_13;
    MY_24 = MY_13;
    MZ_24 = MZ_13;
    MXX_24 = MXX_13;
    MYY_24 = MYY_13;
    MZZ_24 = MZZ_13;
  }
  else {
    compute_multipole_ints(space2, space4,
                           MX_24, MY_24, MZ_24,
                           MXX_24, MYY_24, MZZ_24,
                           MXY_24, MXZ_24, MYZ_24);
    S_24 = compute_overlap_ints(space2, space4);
  }
  if (debug_ >= DefaultPrintThresholds::diagnostics)
    ExEnv::out0() << indent << "Computed overlap and multipole moment integrals" << endl;

  const int nproc = r12world()->world()->msg()->n();
  const int me = r12world()->world()->msg()->me();

  const int n1 = space1->rank();
  const int n2 = space2->rank();
  const int n3 = space3->rank();
  const int n4 = space4->rank();
  const int n12 = n1*n2;
  const int n34 = n3*n4;
  const int n1234 = n12*n34;
  double* r2_array = new double[n1234];
  memset(r2_array,0,n1234*sizeof(double));

  int ij = 0;
  double* ijkl_ptr = r2_array;
  for(int i=0; i<n1; i++)
    for(int j=0; j<n2; j++, ij++) {

    int ij_proc = ij%nproc;
    if (ij_proc != me) {
      ijkl_ptr += n34;
      continue;
    }

    int kl=0;
    for(int k=0; k<n3; k++)
      for(int l=0; l<n4; l++, kl++, ijkl_ptr++) {

        double r2_ik = -1.0*(MXX_13->get_element(i,k) + MYY_13->get_element(i,k) + MZZ_13->get_element(i,k));
        double r2_jl = -1.0*(MXX_24->get_element(j,l) + MYY_24->get_element(j,l) + MZZ_24->get_element(j,l));
        double r11_ijkl = MX_13->get_element(i,k)*MX_24->get_element(j,l) +
          MY_13->get_element(i,k)*MY_24->get_element(j,l) +
          MZ_13->get_element(i,k)*MZ_24->get_element(j,l);
        double S_ik = S_13.get_element(i,k);
        double S_jl = S_24.get_element(j,l);

        double R2_ijkl = r2_ik * S_jl + r2_jl * S_ik - 2.0*r11_ijkl;
        *ijkl_ptr = R2_ijkl;
      }
    }

  r12world()->world()->msg()->sum(r2_array,n1234);

  RefSCDimension dim_ij = new SCDimension(n12);
  RefSCDimension dim_kl = new SCDimension(n34);

  Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();
  RefSCMatrix R2 = local_matrix_kit->matrix(dim_ij, dim_kl);
  R2.assign(r2_array);
  delete[] r2_array;

  // Only task 0 needs R2
  if (me != 0)
    R2.assign(0.0);

  return R2;
}

void
R12IntEval::r2_contrib_to_X_new_()
{
  unsigned int me = r12world()->world()->msg()->me();

  for(int s=0; s<nspincases2(); s++) {

    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
    const Ref<OrbitalSpace>& space1 = GGspace(case1(spincase2));
    const Ref<OrbitalSpace>& space2 = GGspace(case2(spincase2));

    // compute r_{12}^2 operator in act.occ.pair/act.occ.pair basis
    RefSCMatrix R2 = compute_r2_(space1,space2,space1,space2);
    if (spincase2 == AlphaBeta) {
      X_[s].accumulate(R2);
    }
    else {
      // space1 == space2 because it's AlphaAlpha or BetaBeta
      antisymmetrize(X_[s],R2,space1,space1);
    }
  }
}

const Ref<OrbitalSpace>&
R12IntEval::hj_i_P(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return hj_i_P(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = occ_act(spin);
  const Ref<OrbitalSpace>& intspace = r12world()->ribs_space();
  Ref<OrbitalSpace> null;
#if 1 // real H+J
  f_bra_ket(spin,false,true,false,
	    null,
	    hj_i_P_[s],
	    null,
	    extspace,intspace);
#else // use another one-electron matrix instead
  {
    RefSCMatrix X_i_e;
    X_i_e = this->fock(intspace, extspace, AnySpinCase1, 1.0, 0.0, 1.0, 0);
    // r12world()->compute_overlap_ints(intspace, extspace, S_i_e);
    std::string id = extspace->id();  id += "_X(";  id += intspace->id();  id += ")";
    id = ParsedOrbitalSpaceKey::key(id,spin);
    std::string name = "X-weighted space";
    name = prepend_spincase(spin,name);
    hj_i_P_[s] = new OrbitalSpace(id, name, extspace, intspace->coefs()*X_i_e,
                         intspace->basis());
    this->orbital_registry()->add(make_keyspace_pair(hj_i_P_[s]));
  }
#endif
  return hj_i_P_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::hj_i_A(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return hj_i_A(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = occ_act(spin);
  const Ref<OrbitalSpace>& intspace = r12world()->cabs_space(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,true,false,
	    null,
	    hj_i_A_[s],
	    null,
	    extspace,intspace);
  return hj_i_A_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::hj_i_p(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return hj_i_p(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = occ_act(spin);
  const Ref<OrbitalSpace>& intspace = orbs(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,true,false,
	    null,
	    hj_i_p_[s],
	    null,
	    extspace,intspace);
  return hj_i_p_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::hj_i_m(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return hj_i_m(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = occ_act(spin);
  const Ref<OrbitalSpace>& intspace = occ(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,true,false,
	    null,
	    hj_i_m_[s],
	    null,
	    extspace,intspace);
  return hj_i_m_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::hj_i_a(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return hj_i_a(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = occ_act(spin);
  const Ref<OrbitalSpace>& intspace = vir_act(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,true,false,
	    null,
	    hj_i_a_[s],
	    null,
	    extspace,intspace);
  return hj_i_a_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::hj_m_m(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return hj_m_m(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = occ(spin);
  const Ref<OrbitalSpace>& intspace = occ(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,true,false,
        null,
        hj_m_m_[s],
        null,
        extspace,intspace);
  return hj_m_m_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::hj_m_p(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return hj_m_p(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = occ(spin);
  const Ref<OrbitalSpace>& intspace = orbs(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,true,false,
        null,
        hj_m_p_[s],
        null,
        extspace,intspace);
  return hj_m_p_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::hj_a_A(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return hj_a_A(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = vir_act(spin);
  const Ref<OrbitalSpace>& intspace = r12world()->cabs_space(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,true,false,
        null,
        hj_a_A_[s],
        null,
        extspace,intspace);
  return hj_a_A_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::hj_p_P(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return hj_p_P(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = orbs(spin);
  const Ref<OrbitalSpace>& intspace = r12world()->ribs_space();
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,true,false,
	    null,
	    hj_p_P_[s],
	    null,
	    extspace,intspace);
  return hj_p_P_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::hj_p_A(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return hj_p_A(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = orbs(spin);
  const Ref<OrbitalSpace>& intspace = r12world()->cabs_space(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,true,false,
	    null,
	    hj_p_A_[s],
	    null,
	    extspace,intspace);
  return hj_p_A_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::hj_p_p(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return hj_p_p(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = orbs(spin);
  const Ref<OrbitalSpace>& intspace = orbs(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,true,false,
	    null,
	    hj_p_p_[s],
	    null,
	    extspace,intspace);
  return hj_p_p_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::hj_p_m(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return hj_p_m(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = orbs(spin);
  const Ref<OrbitalSpace>& intspace = occ(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,true,false,
	    null,
	    hj_p_m_[s],
	    null,
	    extspace,intspace);
  return hj_p_m_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::hj_p_a(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return hj_p_a(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = orbs(spin);
  const Ref<OrbitalSpace>& intspace = vir_act(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,true,false,
	    null,
	    hj_p_a_[s],
	    null,
	    extspace,intspace);
  return hj_p_a_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::hj_P_P(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return hj_P_P(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = r12world()->ribs_space();
  const Ref<OrbitalSpace>& intspace = r12world()->ribs_space();
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,true,false,
        null,
        hj_P_P_[s],
        null,
        extspace,intspace);
  return hj_P_P_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::K_i_P(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_i_P(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = occ_act(spin);
  const Ref<OrbitalSpace>& intspace = r12world()->ribs_space();
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,false,true,
	    null,
	    null,
	    K_i_P_[s],
	    extspace,intspace);
  return K_i_P_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::K_i_A(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_i_A(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = occ_act(spin);
  const Ref<OrbitalSpace>& intspace = r12world()->cabs_space(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,false,true,
	    null,
	    null,
	    K_i_A_[s],
	    extspace,intspace);
  return K_i_A_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::K_i_p(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_i_p(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = occ_act(spin);
  const Ref<OrbitalSpace>& intspace = orbs(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,false,true,
	    null,
	    null,
	    K_i_p_[s],
	    extspace,intspace);
  return K_i_p_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::K_i_m(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_i_m(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = occ_act(spin);
  const Ref<OrbitalSpace>& intspace = occ(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,false,true,
	    null,
	    null,
	    K_i_m_[s],
	    extspace,intspace);
  return K_i_m_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::K_i_a(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_i_a(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = occ_act(spin);
  const Ref<OrbitalSpace>& intspace = vir_act(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,false,true,
	    null,
	    null,
	    K_i_a_[s],
	    extspace,intspace);
  return K_i_a_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::K_m_a(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_m_a(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = occ(spin);
  const Ref<OrbitalSpace>& intspace = vir_act(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,false,true,
	    null,
	    null,
	    K_m_a_[s],
	    extspace,intspace);
  return K_m_a_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::K_a_a(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_a_a(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = vir_act(spin);
  const Ref<OrbitalSpace>& intspace = vir_act(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,false,true,
	    null,
	    null,
	    K_a_a_[s],
	    extspace,intspace);
  return K_a_a_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::K_a_p(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_a_p(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = vir_act(spin);
  const Ref<OrbitalSpace>& intspace = orbs(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,false,true,
        null,
        null,
        K_a_p_[s],
        extspace,intspace);
  return K_a_p_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::K_a_P(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_a_P(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = vir_act(spin);
  const Ref<OrbitalSpace>& intspace = r12world()->ribs_space();
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,false,true,
        null,
        null,
        K_a_P_[s],
        extspace,intspace);
  return K_a_P_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::K_p_P(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_p_P(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = orbs(spin);
  const Ref<OrbitalSpace>& intspace = r12world()->ribs_space();
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,false,true,
	    null,
	    null,
	    K_p_P_[s],
	    extspace,intspace);
  return K_p_P_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::K_p_A(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_p_A(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = orbs(spin);
  const Ref<OrbitalSpace>& intspace = r12world()->cabs_space(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,false,true,
	    null,
	    null,
	    K_p_A_[s],
	    extspace,intspace);
  return K_p_A_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::K_p_p(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_p_p(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = orbs(spin);
  const Ref<OrbitalSpace>& intspace = orbs(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,false,true,
	    null,
	    null,
	    K_p_p_[s],
	    extspace,intspace);
  return K_p_p_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::K_p_m(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_p_m(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = orbs(spin);
  const Ref<OrbitalSpace>& intspace = occ(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,false,true,
	    null,
	    null,
	    K_p_m_[s],
	    extspace,intspace);
  return K_p_m_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::K_p_a(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_p_a(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = orbs(spin);
  const Ref<OrbitalSpace>& intspace = vir_act(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,false,true,
	    null,
	    null,
	    K_p_a_[s],
	    extspace,intspace);
  return K_p_a_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::K_A_P(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_A_P(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = r12world()->cabs_space(spin);
  const Ref<OrbitalSpace>& intspace = r12world()->ribs_space();
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,false,false,true,
        null,
        null,
        K_A_P_[s],
        extspace,intspace);
  return K_A_P_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::K_P_P(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_P_P(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = r12world()->ribs_space();
  const Ref<OrbitalSpace>& intspace = r12world()->ribs_space();

#define TEST_FOCKBUILD_RUNTIME 0
#if TEST_FOCKBUILD_RUNTIME
  {
    Ref<FockBuildRuntime> fb_rtime = fockbuild_runtime();
    const SpinCase1 realspin = r12world()->refwfn()->refwfn()->spin_polarized() ? spin : AnySpinCase1;

    // get AO basis matrix
    const Ref<OrbitalSpace>& ribs = intspace;
    Ref<OrbitalSpace> ribs_ao = this->ao_registry()->value(ribs->basis());
    Ref<OrbitalSpace> null;
    Ref<OrbitalSpace> K1;
    f_bra_ket(spin,false,false,true,
          null,
          null,
          K1,
          ribs_ao,ribs_ao);
    K1->coefs().print("K in CABS+/CABS+ AO basis using old transform-based build");

    const std::string kkey_ao = ParsedOneBodyIntKey::key(ribs_ao->id(),ribs_ao->id(),std::string("K"),realspin);
    RefSCMatrix K2 = fb_rtime->get(kkey_ao);
    K2.print("K in CABS+/CABS+ AO basis using new build");

    const std::string kkey = ParsedOneBodyIntKey::key(ribs->id(),ribs->id(),std::string("K"),realspin);
    RefSCMatrix K3 = fb_rtime->get(kkey);
    K3.print("K in CABS+/CABS+ basis using new build");

    Ref<OrbitalSpace> orbs = orbs(spin);
    Ref<OrbitalSpace> obs_ao = this->ao_registry()->value(orbs->basis());
    Ref<OrbitalSpace> K4;
    f_bra_ket(spin,false,false,true,
          null,
          null,
          K4,
          obs_ao,obs_ao);
    K4->coefs().print("K in OBS/OBS AO basis using old transform-based build");

    {
    const std::string kkey_ao = ParsedOneBodyIntKey::key(obs_ao->id(),obs_ao->id(),std::string("K"),realspin);
    RefSCMatrix K5 = fb_rtime->get(kkey_ao);
    K5.print("K in OBS/OBS AO basis using new build");
    }

  }
#endif

  Ref<OrbitalSpace> null;
  f_bra_ket(spin,true,false,true,
	    F_P_P_[s],
	    null,
	    K_P_P_[s],
	    extspace,intspace);
  return K_P_P_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::F_P_P(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_P_P(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = r12world()->ribs_space();
  const Ref<OrbitalSpace>& intspace = r12world()->ribs_space();
  Ref<OrbitalSpace> null;
#define TEST_FOCK_BUILDER 0
  f_bra_ket(spin,true,false,
#if !TEST_FOCK_BUILDER
            true,
#else
            false,
#endif
	    F_P_P_[s],
	    null,
#if !TEST_FOCK_BUILDER
	    K_P_P_[s],
#else
	    null,
#endif
	    extspace,intspace);
  return F_P_P_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::h_P_P(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return h_P_P(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  if (h_P_P_[s]) return h_P_P_[s];

  const Ref<OrbitalSpace>& extspace = r12world()->ribs_space();
  const Ref<OrbitalSpace>& intspace = r12world()->ribs_space();

  RefSCMatrix h_i_e = fock(intspace, extspace, spin, 0.0, 0.0, 1.0);
  if (debug_ >= DefaultPrintThresholds::allN2) {
    std::string label("h matrix in ");
    label += intspace->id();
    label += "/";
    label += extspace->id();
    label += " basis";
    h_i_e.print(label.c_str());
  }
  std::string id = extspace->id();
  id += "_h(";
  id += intspace->id();
  id += ")";
  id = ParsedOrbitalSpaceKey::key(id, spin);
  std::string name = "h-weighted space";
  name = prepend_spincase(spin, name);

  h_P_P_[s] = new OrbitalSpace(id, name, extspace, intspace->coefs() * h_i_e,
                               intspace->basis());

  const Ref<OrbitalSpaceRegistry>& idxreg = this->orbital_registry();
  idxreg->add(make_keyspace_pair(h_P_P_[s]));

  return h_P_P_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::gamma_p_p(SpinCase1 S) {
  MPQC_ASSERT(!r12world()->spinadapted());
  if (!spin_polarized() && S == Beta)
    return gamma_p_p(Alpha);

  if (gamma_p_p_[S].null()) {
    const Ref<OrbitalSpace>& extspace = this->orbs(S);
    const Ref<OrbitalSpace>& intspace = this->orbs(S);
    std::string id = extspace->id();
    id += "_gamma";
    if (spin_polarized())
      id += (S == Alpha ? "[A]" : "[B]");
    id += "(";  id += intspace->id();  id += ")";
    std::string name = "gamma-weighted space";
    gamma_p_p_[S] = new OrbitalSpace(id, name, extspace, intspace->coefs() * this->ordm(S),
                                     intspace->basis());
    this->orbital_registry()->add(make_keyspace_pair(gamma_p_p_[S]));
  }
  return gamma_p_p_[S];
}



const Ref<OrbitalSpace>&
R12IntEval::gamma_p_p_av() {
  MPQC_ASSERT(r12world()->spinadapted());
  if (gamma_p_p_[Alpha].null()) {
    const Ref<OrbitalSpace>& extspace = this->orbs(Alpha);
    const Ref<OrbitalSpace>& intspace = this->orbs(Alpha);
    std::string id = extspace->id();  id += "_gamma(";  id += intspace->id();  id += ")";
    std::string name = "gamma-weighted space";
    gamma_p_p_[Alpha] = new OrbitalSpace(id, name, extspace, intspace->coefs() * this->ordm_av(),
                                     intspace->basis());
    this->orbital_registry()->add(make_keyspace_pair(gamma_p_p_[Alpha]));
  }
  return gamma_p_p_[Alpha];
}

const Ref<OrbitalSpace>&
R12IntEval::gamma_m_m_av() {
  MPQC_ASSERT(r12world()->spinadapted());
  if (gamma_m_m_[Alpha].null()) {
    const Ref<OrbitalSpace>& extspace = this->occ(Alpha);
    const Ref<OrbitalSpace>& intspace = this->occ(Alpha);
    std::string id = extspace->id();  id += "_gamma(";  id += intspace->id();  id += ")";
    std::string name = "gamma-weighted space";
    gamma_m_m_[Alpha] = new OrbitalSpace(id, name, extspace, intspace->coefs() * this->ordm_occ_av(),
                                     intspace->basis());
    this->orbital_registry()->add(make_keyspace_pair(gamma_m_m_[Alpha]));
  }
  return gamma_m_m_[Alpha];
}

const Ref<OrbitalSpace>&
R12IntEval::gammaFgamma_p_p(SpinCase1 S) {
  MPQC_ASSERT(!r12world()->spinadapted());
  if (!spin_polarized() && S == Beta)
    return gammaFgamma_p_p(Alpha);

  if (gammaFgamma_p_p_[S].null()) {
    const Ref<OrbitalSpace>& extspace = this->orbs(S);
    const Ref<OrbitalSpace>& intspace = this->orbs(S);
    RefSCMatrix F_i_e = fock(intspace,extspace,S,1.0,1.0);
    std::string id = extspace->id();  id += "_gFg";
    if (spin_polarized())
      id += (S == Alpha ? "[A]" : "[B]");
    id += "(";  id += intspace->id();  id += ")";
    std::string name = "gammaFgamma-weighted space";
    gammaFgamma_p_p_[S] = new OrbitalSpace(id, name, extspace, intspace->coefs() * this->ordm(S) * F_i_e * this->ordm(S),
                                           intspace->basis());
    this->orbital_registry()->add(make_keyspace_pair(gammaFgamma_p_p_[S]));
  }
  return gammaFgamma_p_p_[S];
}

const Ref<OrbitalSpace>&
R12IntEval::gammaFgamma_p_p() {
  MPQC_ASSERT(r12world()->spinadapted());
  if (gammaFgamma_p_p_[Alpha].null()) {
    const Ref<OrbitalSpace>& extspace = this->orbs(Alpha);
    const Ref<OrbitalSpace>& intspace = this->orbs(Alpha);
    RefSCMatrix F_i_e = fock(intspace,extspace,Alpha,1.0,1.0);
    std::string id = extspace->id();  id += "_gFg";
    id += "(";  id += intspace->id();  id += ")";
    std::string name = "gammaFgamma-weighted space";
    gammaFgamma_p_p_[Alpha] = new OrbitalSpace(id, name, extspace, intspace->coefs() * this->ordm_av() * F_i_e * this->ordm_av(),
                                           intspace->basis());
    this->orbital_registry()->add(make_keyspace_pair(gammaFgamma_p_p_[Alpha]));
  }
  return gammaFgamma_p_p_[Alpha];
}

const Ref<OrbitalSpace>&
R12IntEval::Fgamma_p_P(SpinCase1 S) {
  MPQC_ASSERT(!r12world()->spinadapted());
  if (!spin_polarized() && S == Beta)
    return(Fgamma_p_P(Alpha));
  if(Fgamma_p_P_[S].null()) {
    const Ref<OrbitalSpace>& extspace = this->orbs(S);
    const Ref<OrbitalSpace>& intspace = r12world()->ribs_space();
    RefSCMatrix F_i_e = fock(intspace,extspace,S,1.0,1.0);
    std::string id = extspace->id();  id += "_Fg";
    if (spin_polarized())
      id += (S == Alpha ? "[A]" : "[B]");
    id += "(";  id += intspace->id();  id += ")";
    std::string name = "Fgamma-weighted space";
    Fgamma_p_P_[S] = new OrbitalSpace(id, name, extspace, intspace->coefs() * F_i_e * this->ordm(S),
                                      intspace->basis());
    this->orbital_registry()->add(make_keyspace_pair(Fgamma_p_P_[S]));
  }
  return Fgamma_p_P_[S];
}


const Ref<OrbitalSpace>&
R12IntEval::Fgamma_p_P() {
  MPQC_ASSERT(r12world()->spinadapted());
  if(Fgamma_p_P_[Alpha].null()) {
    const Ref<OrbitalSpace>& extspace = this->orbs(Alpha);
    const Ref<OrbitalSpace>& intspace = r12world()->ribs_space();
    RefSCMatrix F_i_e = fock(intspace,extspace,Alpha,1.0,1.0);
    std::string id = extspace->id();  id += "_Fg(";  id += intspace->id();  id += ")";
    std::string name = "Fgamma-weighted space";
    Fgamma_p_P_[Alpha] = new OrbitalSpace(id, name, extspace, intspace->coefs() * F_i_e * this->ordm_av(),
                                      intspace->basis());
    this->orbital_registry()->add(make_keyspace_pair(Fgamma_p_P_[Alpha]));
  }
  return Fgamma_p_P_[Alpha];
}



Ref<OrbitalSpace> R12IntEval::obtensor_p_A(const RefSCMatrix &obtensor,SpinCase1 S) {
  const Ref<OrbitalSpace>& extspace = this->orbs(S);
  const Ref<OrbitalSpace>& intspace = r12world()->cabs_space(S);

  RefSCDimension dimA = intspace->coefs()->coldim();
  RefSCDimension dimp = extspace->coefs()->coldim();
  RefSCMatrix obtensor_blkd = intspace->coefs()->kit()->matrix(dimA,dimp);
  for(int i=0; i<dimA.n(); i++) {
    for(int j=0; j<dimp.n(); j++) {
      obtensor_blkd.set_element(i,j,obtensor.get_element(i,j));
    }
  }

  std::string id = extspace->id();  id += "_obt(";  id += intspace->id();  id += ")";
  std::string name = "obtensor-weighted space";
  Ref<OrbitalSpace> obtensor_moidx = new OrbitalSpace(id, name, extspace, intspace->coefs()*obtensor_blkd,
                                                      intspace->basis());
  return(obtensor_moidx);
}

const Ref<OrbitalSpace>&
R12IntEval::F_p_A(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_p_A(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = orbs(spin);
  const Ref<OrbitalSpace>& intspace = r12world()->cabs_space(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,true,false,false,
	    F_p_A_[s],
	    null,
	    null,
	    extspace,intspace);
  return F_p_A_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::F_p_p(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_p_p(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = orbs(spin);
  const Ref<OrbitalSpace>& intspace = orbs(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,true,false,false,
	    F_p_p_[s],
	    null,
	    null,
	    extspace,intspace);
  return F_p_p_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::F_p_m(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_p_m(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = orbs(spin);
  const Ref<OrbitalSpace>& intspace = occ(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,true,false,false,
        F_p_m_[s],
        null,
        null,
        extspace,intspace);
  return F_p_m_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::F_p_a(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_p_a(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = orbs(spin);
  const Ref<OrbitalSpace>& intspace = vir_act(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,true,false,false,
        F_p_a_[s],
        null,
        null,
        extspace,intspace);
  return F_p_a_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::F_m_m(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_m_m(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = occ(spin);
  const Ref<OrbitalSpace>& intspace = occ(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,true,false,false,
	    F_m_m_[s],
	    null,
	    null,
	    extspace,intspace);
  return F_m_m_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::F_m_a(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_m_a(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = occ(spin);
  const Ref<OrbitalSpace>& intspace = vir_act(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,true,false,false,
	    F_m_a_[s],
	    null,
	    null,
	    extspace,intspace);
  return F_m_a_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::F_m_p(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_m_p(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = occ(spin);
  const Ref<OrbitalSpace>& intspace = orbs(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,true,false,false,
        F_m_p_[s],
        null,
        null,
        extspace,intspace);
  return F_m_p_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::F_m_P(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_m_P(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = occ(spin);
  const Ref<OrbitalSpace>& intspace = r12world()->ribs_space();
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,true,false,false,
	    F_m_P_[s],
	    null,
	    null,
	    extspace,intspace);
  return F_m_P_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::F_gg_P(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_gg_P(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = ggspace(spin);
  const Ref<OrbitalSpace>& intspace = r12world()->ribs_space();
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,true,false,false,
        F_gg_P_[s],
        null,
        null,
        extspace,intspace);
  return F_gg_P_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::F_m_A(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_m_A(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = occ(spin);
  const Ref<OrbitalSpace>& intspace = r12world()->cabs_space(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,true,false,false,
	    F_m_A_[s],
	    null,
	    null,
	    extspace,intspace);
  return F_m_A_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::F_i_A(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_i_A(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = occ_act(spin);
  const Ref<OrbitalSpace>& intspace = r12world()->cabs_space(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,true,false,false,
	    F_i_A_[s],
	    null,
	    null,
	    extspace,intspace);
  return F_i_A_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::F_i_p(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_i_p(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = occ_act(spin);
  const Ref<OrbitalSpace>& intspace = orbs(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,true,false,false,
        F_i_p_[s],
        null,
        null,
        extspace,intspace);
  return F_i_p_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::F_i_m(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_i_m(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = occ_act(spin);
  const Ref<OrbitalSpace>& intspace = occ(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,true,false,false,
        F_i_m_[s],
        null,
        null,
        extspace,intspace);
  return F_i_m_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::F_i_a(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_i_a(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = occ_act(spin);
  const Ref<OrbitalSpace>& intspace = vir_act(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,true,false,false,
        F_i_a_[s],
        null,
        null,
        extspace,intspace);
  return F_i_a_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::F_i_P(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_i_P(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = occ_act(spin);
  const Ref<OrbitalSpace>& intspace = r12world()->ribs_space();
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,true,false,false,
        F_i_P_[s],
        null,
        null,
        extspace,intspace);
  return F_i_P_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::F_a_a(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_a_a(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = vir_act(spin);
  const Ref<OrbitalSpace>& intspace = vir_act(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,true,false,false,
	    F_a_a_[s],
	    null,
	    null,
	    extspace,intspace);
  return F_a_a_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::F_a_A(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_a_A(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = vir_act(spin);
  const Ref<OrbitalSpace>& intspace = r12world()->cabs_space(spin);
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,true,false,false,
	    F_a_A_[s],
	    null,
	    null,
	    extspace,intspace);
  return F_a_A_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::J_i_p(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return J_i_p(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  if (J_i_p_[s]) return J_i_p_[s];

  const Ref<OrbitalSpace>& extspace = occ_act(spin);
  const Ref<OrbitalSpace>& intspace = orbs(spin);

  RefSCMatrix J_i_e = fock(intspace, extspace, spin, 1.0, 0.0, 0.0);
  if (debug_ >= DefaultPrintThresholds::allN2) {
    std::string label("J matrix in ");
    label += intspace->id();
    label += "/";
    label += extspace->id();
    label += " basis";
    J_i_e.print(label.c_str());
  }
  std::string id = extspace->id();
  id += "_J(";
  id += intspace->id();
  id += ")";
  id = ParsedOrbitalSpaceKey::key(id, spin);
  std::string name = "J-weighted space";
  name = prepend_spincase(spin, name);

  J_i_p_[s] = new OrbitalSpace(id, name, extspace, intspace->coefs() * J_i_e,
                                intspace->basis());

  const Ref<OrbitalSpaceRegistry>& idxreg = this->orbital_registry();
  idxreg->add(make_keyspace_pair(J_i_p_[s]));

  return J_i_p_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::J_i_P(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return J_i_P(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  if (J_i_P_[s]) return J_i_P_[s];

  const Ref<OrbitalSpace>& extspace = occ_act(spin);
  const Ref<OrbitalSpace>& intspace = r12world()->ribs_space();

  RefSCMatrix J_i_e = fock(intspace, extspace, spin, 1.0, 0.0, 0.0);
  if (debug_ >= DefaultPrintThresholds::allN2) {
    std::string label("J matrix in ");
    label += intspace->id();
    label += "/";
    label += extspace->id();
    label += " basis";
    J_i_e.print(label.c_str());
  }
  std::string id = extspace->id();
  id += "_J(";
  id += intspace->id();
  id += ")";
  id = ParsedOrbitalSpaceKey::key(id, spin);
  std::string name = "J-weighted space";
  name = prepend_spincase(spin, name);

  J_i_P_[s] = new OrbitalSpace(id, name, extspace, intspace->coefs() * J_i_e,
                                intspace->basis());

  const Ref<OrbitalSpaceRegistry>& idxreg = this->orbital_registry();
  idxreg->add(make_keyspace_pair(J_i_P_[s]));

  return J_i_P_[s];
}

const Ref<OrbitalSpace>&
R12IntEval::J_P_P(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return J_P_P(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  if (J_P_P_[s]) return J_P_P_[s];

  const Ref<OrbitalSpace>& extspace = r12world()->ribs_space();
  const Ref<OrbitalSpace>& intspace = r12world()->ribs_space();

  RefSCMatrix J_i_e = fock(intspace, extspace, spin, 1.0, 0.0, 0.0);
  if (debug_ >= DefaultPrintThresholds::allN2) {
    std::string label("J matrix in ");
    label += intspace->id();
    label += "/";
    label += extspace->id();
    label += " basis";
    J_i_e.print(label.c_str());
  }
  std::string id = extspace->id();
  id += "_J(";
  id += intspace->id();
  id += ")";
  id = ParsedOrbitalSpaceKey::key(id, spin);
  std::string name = "J-weighted space";
  name = prepend_spincase(spin, name);

  J_P_P_[s] = new OrbitalSpace(id, name, extspace, intspace->coefs() * J_i_e,
                                intspace->basis());

  const Ref<OrbitalSpaceRegistry>& idxreg = this->orbital_registry();
  idxreg->add(make_keyspace_pair(J_P_P_[s]));

  return J_P_P_[s];
}

const Ref<OrbitalSpace>& R12IntEval::F_A_A(SpinCase1 spin) {
  if (!spin_polarized() && spin == Beta)
    return F_A_A(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = r12world()->cabs_space(spin);
  const Ref<OrbitalSpace>& intspace = r12world()->cabs_space(spin);

  Ref<OrbitalSpace> null;
  f_bra_ket(spin,true,false,false,
        F_A_A_[s],
        null,
        null,
        extspace,intspace);
  return(F_A_A_[s]);
}

const Ref<OrbitalSpace>& R12IntEval::F_p_P(SpinCase1 spin) {
  if (!spin_polarized() && spin == Beta)
    return(F_p_P(Alpha));

  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<OrbitalSpace>& extspace = orbs(spin);
  const Ref<OrbitalSpace>& intspace = r12world()->ribs_space();
  Ref<OrbitalSpace> null;
  f_bra_ket(spin,true,false,false,
        F_p_P_[s],
        null,
        null,
        extspace,intspace);
  return(F_p_P_[s]);
}

void
R12IntEval::f_bra_ket(
    SpinCase1 spin,
    bool make_F,
    bool make_hJ,
    bool make_K,
    Ref<OrbitalSpace>& F,
    Ref<OrbitalSpace>& hJ,
    Ref<OrbitalSpace>& K,
    const Ref<OrbitalSpace>& extspace,
    const Ref<OrbitalSpace>& intspace
    )
{
  const Ref<OrbitalSpaceRegistry>& idxreg = this->orbital_registry();

  const unsigned int s = static_cast<unsigned int> (spin);
  const bool not_yet_computed = (make_F && F.null()) || (make_hJ && hJ.null())
      || (make_K && K.null());
  if (not_yet_computed) {

    const int dk = this->dk();

#if 0
    ExEnv::out0() << indent << "make_F   = " << (make_F ? "true" : "false") << endl;
    ExEnv::out0() << indent << "make_hJ  = " << (make_hJ ? "true" : "false") << endl;
    ExEnv::out0() << indent << "make_K   = " << (make_K ? "true" : "false") << endl;
#endif

    RefSCMatrix hJ_i_e;
    if (make_hJ && hJ.null()) {
      hJ_i_e = fock(intspace, extspace, spin, 1.0, 0.0);
      if (debug_ >= DefaultPrintThresholds::allN2) {
        std::string label("(h+J) matrix in ");
        label += intspace->id();
        label += "/";
        label += extspace->id();
        label += " basis";
        hJ_i_e.print(label.c_str());
      }

      std::string id = extspace->id();
      id += "_hJ(";
      id += intspace->id();
      id += ")";
      id = ParsedOrbitalSpaceKey::key(id, spin);
      std::string name = "(h+J)-weighted space";
      name = prepend_spincase(spin, name);

      hJ = new OrbitalSpace(id, name, extspace, intspace->coefs() * hJ_i_e,
                            intspace->basis());
      idxreg->add(make_keyspace_pair(hJ));
    }

    RefSCMatrix K_i_e;
    if (make_K && K.null()) {
      if (!USE_FOCKBUILD)
        K_i_e = exchange_(spin, intspace, extspace);
      else {
        K_i_e = fock(intspace, extspace, spin, 0.0, 1.0, 0.0);
        K_i_e.scale(-1.0);
      }
      if (debug_ >= DefaultPrintThresholds::allN2) {
        std::string label;
        label += "K matrix in ";
        label += intspace->id();
        label += "/";
        label += extspace->id();
        label += " basis";
        K_i_e.print(label.c_str());
      }

      std::string id = extspace->id();
      id += "_K(";
      id += intspace->id();
      id += ")";
      id = ParsedOrbitalSpaceKey::key(id, spin);
      std::string name = "K-weighted space";
      name = prepend_spincase(spin, name);
      K = new OrbitalSpace(id, name, extspace, intspace->coefs() * K_i_e,
                           intspace->basis());
      idxreg->add(make_keyspace_pair(K));
    }

    if (make_F && F.null()) {
      RefSCMatrix F_i_e;
      if (make_hJ) {
        if (make_K) {
          F_i_e = K_i_e.clone();
          F_i_e.assign(K_i_e);
          F_i_e.scale(-1.0);
          F_i_e.accumulate(hJ_i_e);
        } else {
          if (!USE_FOCKBUILD) {
            F_i_e = exchange_(spin, intspace, extspace);
            F_i_e.scale(-1.0);
          } else {
            F_i_e = fock(intspace, extspace, spin, 0.0, 1.0, 0.0);
          }
          F_i_e.accumulate(hJ_i_e);
        }
      } else {
        if (make_K) {
          F_i_e = K_i_e.clone();
          F_i_e.assign(K_i_e);
          F_i_e.scale(-1.0);
          RefSCMatrix hJ_i_e = fock(intspace, extspace, spin, 1.0, 0.0);
          F_i_e.accumulate(hJ_i_e);
        } else {
          F_i_e = fock(intspace, extspace, spin, 1.0, 1.0);
        }
      }
      if (debug_ >= DefaultPrintThresholds::allN2) {
        std::string label("F matrix in ");
        label += intspace->id();
        label += "/";
        label += extspace->id();
        label += " basis";
        F_i_e.print(label.c_str());
      }

      std::string id = extspace->id();
      id += "_F(";
      id += intspace->id();
      id += ")";
      id = ParsedOrbitalSpaceKey::key(id, spin);
      std::string name = "F-weighted space";
      name = prepend_spincase(spin, name);
      F = new OrbitalSpace(id, name, extspace, intspace->coefs() * F_i_e,
                           intspace->basis());
      idxreg->add(make_keyspace_pair(F));
    }
  }
}

void
R12IntEval::compute()
{
  if (evaluated_)
    return;

#if 0
  {
    ExEnv::out0() << indent << "Trying out MPQC3-based R12 code" << std::endl;
    //V_diag_ta();
    gf2_r12();
  }
#endif

  init_intermeds_();

  // different expressions hence codepaths depending on relationship between OBS, VBS, and RIBS
  // compare these basis sets here
  const bool obs_eq_vbs = r12world()->obs_eq_vbs();
  const bool obs_eq_ribs = r12world()->obs_eq_ribs();
  const R12Technology::ABSMethod absmethod = r12world()->r12tech()->abs_method();
  // is CABS space empty?
  const bool cabs_empty = obs_eq_vbs && obs_eq_ribs;
  const bool vir_empty = vir(Alpha)->rank()==0 || vir(Beta)->rank()==0;

  Ref<R12Technology::NullCorrelationFactor> nocorrptr; nocorrptr << corrfactor();
  // if explicit correlation -- compute F12 theory intermediates
  if (nocorrptr.null()) {

    if (debug_ >= DefaultPrintThresholds::O4) {
      globally_sum_intermeds_();
      for(int s=0; s<nspincases2(); s++) {
        V_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"V(diag) contribution").c_str());
        X_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"X(diag) contribution").c_str());
        B_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"B(diag) contribution").c_str());
      }
    }

    if (obs_eq_vbs) {
      if(r12world()->r12tech()->ansatz()->projector()==R12Technology::Projector_1) {
        if (!r12world()->sdref())
          throw FeatureNotImplemented("Projector 1 cannot be used with multi-configuration references",
                                      __FILE__, __LINE__);
        R12IntEval::contrib_to_VXB_c_ansatz1_();
      }
      else {
        if (r12world()->spinadapted()) // to compute in spin-adapted fashion assume general reference
                                       // would be more efficient to have specialized single-determinant spin-adapted version, but not there yet
          contrib_to_VX_GenRefansatz2_spinfree_();
        else if (r12world()->sdref()) // single-determinant reference
          contrib_to_VXB_a_();
        else // spin-orbital general-reference case
          contrib_to_VX_GenRefansatz2_();
      }
    }
    else {
      contrib_to_VXB_a_vbsneqobs_();
    }

    // In stdapprox A', A'', and B add single-commutator contributions due to the relativistic terms
    // in stdapprox C these are computed via RI
    if (this->dk() > 0 && stdapprox() != R12Technology::StdApprox_C) {
      contrib_to_B_DKH_a_();
    }

    // Contribution from X to B in approximation A'' is more complicated than in other methods
    // because exchange is completely skipped
    if (stdapprox() == R12Technology::StdApprox_App) {
      compute_BApp_();
    }
    // whereas other methods that include X (A' and B) can use the simple fX form
    else {
      if (stdapprox() == R12Technology::StdApprox_Ap ||
          stdapprox() == R12Technology::StdApprox_B) {
        compute_B_fX_();
      }
    }

    // This is app B contribution to B -- only valid for projector 2
    if ((stdapprox() == R12Technology::StdApprox_B) && ansatz()->projector() ==
         R12Technology::Projector_2 && !r12world()->r12tech()->omit_B()) {
      compute_BB_();
      if (debug_ >= DefaultPrintThresholds::O4)
        for (int s=0; s<nspincases2(); s++)
          BB_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"B(app. B) contribution").c_str());
    }

    if (!r12world()->r12tech()->omit_B()) {
      if (stdapprox() == R12Technology::StdApprox_C) {
        if(r12world()->r12tech()->ansatz()->projector()==R12Technology::Projector_1){
          compute_BC_ansatz1_();
          if(!r12world()->sdref())
            throw FeatureNotImplemented("Projector 1 cannot be used with multi-configuration references",
                                        __FILE__,__LINE__);
        }
        else {
          if(r12world_->spinadapted()) // see comments above -- treat all spinadapted computations as general-reference
            compute_BC_GenRefansatz2_spinfree();
          else if(r12world()->sdref()) // spin-orbital single-determinant reference
            compute_BC_();
          else // spin-orbital general reference
            compute_BC_GenRefansatz2_();
        }
        if (debug_ >= DefaultPrintThresholds::O4)
          for(int s=0; s<nspincases2(); s++)
            B_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"B(app. C) intermediate").c_str());
      }
      if (stdapprox() == R12Technology::StdApprox_Cp) {
        compute_BCp_();
        if (debug_ >= DefaultPrintThresholds::O4)
          for(int s=0; s<nspincases2(); s++)
            B_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"B(app. C') intermediate").c_str());
      }
    }

#if INCLUDE_COUPLING_CODE
    const bool nonzero_coupling_terms = coupling() && !cabs_empty && !vir_empty;
    if (nonzero_coupling_terms) {
      compute_A();
      if (debug_ >= DefaultPrintThresholds::O2N2) {
        for(int s=0; s<nspincases2(); s++) {
          const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
          std::string label = prepend_spincase(spincase2,"T2 matrix");
          amps()->T2(spincase2).print(label.c_str());
          label = prepend_spincase(spincase2,"F12(vv) matrix");
          amps()->Fvv(spincase2).print(label.c_str());
          label = prepend_spincase(spincase2,"A matrix");
          A_[s].print(label.c_str());
        }
      }

      AT2_contrib_to_V_();
    }
#endif

#if INCLUDE_EBC_CODE
    const bool nonzero_ebc_terms = !ebc() && !cabs_empty && !vir_empty && !r12world()->r12tech()->omit_B();
    if (nonzero_ebc_terms) {
      // EBC contribution to B only appears in commutator-based projector 2 case
      if (ansatz()->projector() == R12Technology::Projector_2 &&
          (stdapprox() == R12Technology::StdApprox_Ap ||
           stdapprox() == R12Technology::StdApprox_App ||
           stdapprox() == R12Technology::StdApprox_B)
          ) {
        AF12_contrib_to_B_();
      }
    }
#endif

#if INCLUDE_GBC_CODE
    // GBC contribution to B only appears in non-StdApproxC projector 2 case
    const bool nonzero_gbc_terms = !gbc() && !cabs_empty &&
                                   ansatz()->projector() == R12Technology::Projector_2 &&
                                   stdapprox() != R12Technology::StdApprox_C &&
                                   !r12world()->r12tech()->omit_B();
    if (nonzero_gbc_terms) {
      // These functions assume that virtuals are expanded in the same basis
      // as the occupied orbitals
      if (!obs_eq_vbs)
        throw std::runtime_error("R12IntEval::compute() -- gbc=false is only supported when basis_vir == basis");

      compute_B_gbc_();
    }
#endif

  } // nontrivial correlation factor

  // Distribute the final intermediates to every node
  globally_sum_intermeds_(true);

  evaluated_ = true;
  return;
}

void
R12IntEval::globally_sum_scmatrix_(RefSCMatrix& A, bool to_all_tasks, bool to_average)
{
  Ref<MessageGrp> msg = r12world()->world()->msg();
  unsigned int ntasks = msg->n();
  // If there's only one task then there's nothing to do
  if (ntasks == 1)
    return;

  const int nelem = A.ncol() * A.nrow();
  double *A_array = new double[nelem];
  A.convert(A_array);
  if (to_all_tasks)
    msg->sum(A_array,nelem,0,-1);
  else
    msg->sum(A_array,nelem,0,0);
  A.assign(A_array);
  if (to_average)
    A.scale(1.0/(double)ntasks);
  if (!to_all_tasks && msg->me() != 0)
    A.assign(0.0);

  delete[] A_array;
}

void
R12IntEval::globally_sum_scvector_(RefSCVector& A, bool to_all_tasks, bool to_average)
{
  Ref<MessageGrp> msg = r12world()->world()->msg();
  unsigned int ntasks = msg->n();
  // If there's only one task then there's nothing to do
  if (ntasks == 1)
    return;

  const int nelem = A.dim().n();
  double *A_array = new double[nelem];
  A.convert(A_array);
  if (to_all_tasks)
    msg->sum(A_array,nelem,0,-1);
  else
    msg->sum(A_array,nelem,0,0);
  A.assign(A_array);
  if (to_average)
    A.scale(1.0/(double)ntasks);
  if (!to_all_tasks && msg->me() != 0)
    A.assign(0.0);

  delete[] A_array;
}

void
R12IntEval::globally_sum_intermeds_(bool to_all_tasks)
{
  for(int s=0; s<nspincases2(); s++) {
    globally_sum_scmatrix_(V_[s],to_all_tasks);
    globally_sum_scmatrix_(X_[s],to_all_tasks);
    globally_sum_scmatrix_(B_[s],to_all_tasks);
    if (BB_[s]) globally_sum_scmatrix_(BB_[s],to_all_tasks);
    if (A_[s]) globally_sum_scmatrix_(A_[s],to_all_tasks);
  }

  if (debug_ >= DefaultPrintThresholds::diagnostics) {
    ExEnv::out0() << indent << "Collected contributions to the intermediates from all tasks";
    if (to_all_tasks)
      ExEnv::out0() << " and distributed to every task" << endl;
    else
      ExEnv::out0() << " on task 0" << endl;
  }
}

const Ref<OrbitalSpace>&
R12IntEval::occ_act(SpinCase1 S) const
{
  return r12world()->refwfn()->occ_act_sb(S);
}

const Ref<OrbitalSpace>&
R12IntEval::occ(SpinCase1 S) const
{
  return r12world()->refwfn()->occ_sb(S);
}

const Ref<OrbitalSpace>&
R12IntEval::vir_act(SpinCase1 S) const
{
  return r12world()->refwfn()->uocc_act_sb(S);
}

const Ref<OrbitalSpace>&
R12IntEval::vir(SpinCase1 S) const
{
  return r12world()->refwfn()->uocc_sb(S);
}

const Ref<OrbitalSpace>&
R12IntEval::orbs(SpinCase1 S) const
{
  return r12world()->refwfn()->orbs_sb(S);
}


const Ref<OrbitalSpace>&
R12IntEval::GGspace(SpinCase1 S) const {
  switch(r12world()->r12tech()->ansatz()->orbital_product_GG()) {
  case R12Technology::OrbProdGG_ij:
  return(occ_act(S));
  case R12Technology::OrbProdGG_pq:
  return(orbs(S));
  default:
  throw ProgrammingError("R12IntEval::GGspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
  }
}

const Ref<OrbitalSpace>& R12IntEval::ggspace(SpinCase1 S) const {
  switch(r12world()->r12tech()->ansatz()->orbital_product_gg()) {
  case R12Technology::OrbProdgg_ij:
  return(occ_act(S));
  case R12Technology::OrbProdgg_pq:
  return(orbs(S));
  default:
  throw ProgrammingError("R12IntEval::ggspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
  }
}

RefSymmSCMatrix
R12IntEval::ordm(SpinCase1 S) const {
  if (r12world()->refwfn()->spin_polarized() == false && S == Beta)
    return ordm(Alpha);
  if (ordm_[S])
    return ordm_[S];

  ordm_[S] = r12world()->refwfn()->ordm_orbs_sb(S);
  return ordm_[S];
}

RefSymmSCMatrix
R12IntEval::ordm() const {
  MPQC_ASSERT(r12world()->spinadapted());
  if (ordm_[Alpha])
    return ordm_[Alpha];

  ordm_[Alpha] = r12world()->refwfn()->ordm_orbs_sb(Alpha);
  ordm_[Alpha].accumulate(r12world()->refwfn()->ordm_orbs_sb(Beta));
  return ordm_[Alpha];
}



RefSymmSCMatrix
R12IntEval::ordm_av() const {
  MPQC_ASSERT(r12world()->spinadapted());
  RefSymmSCMatrix ordm_av = r12world()->refwfn()->ordm_orbs_sb(Alpha).copy();
  ordm_av.accumulate(r12world()->refwfn()->ordm_orbs_sb(Beta));
  ordm_av.scale(0.5);
  return ordm_av;
//  RefSymmSCMatrix adm = r12world()->refwfn()->ordm_orbs_sb(Alpha);
//  adm.accumulate(r12world()->refwfn()->ordm_orbs_sb(Beta));
//  adm.scale(0.5);
//  ordm_[Alpha] = adm;
//  return ordm_[Alpha] ;
}

RefSymmSCMatrix
R12IntEval::ordm_occ_av() const {
  MPQC_ASSERT(r12world()->spinadapted());
  RefSymmSCMatrix ordm_av = r12world()->refwfn()->ordm_occ_sb(Alpha).copy();
  ordm_av.accumulate(r12world()->refwfn()->ordm_occ_sb(Beta));
  ordm_av.scale(0.5);
  return ordm_av;
}

bool
R12IntEval::bc() const {
  Ref<SD_RefWavefunction> sdptr; sdptr << r12world()->refwfn();
  if (sdptr) {
    if (!sdptr->obwfn()->spin_unrestricted() && sdptr->spin_polarized())
      return false; // ROHF
    else {
      Ref<GaussianBasisSet> obs = sdptr->occ_sb(Alpha)->basis();
      Ref<GaussianBasisSet> vbs = sdptr->uocc_sb(Alpha)->basis();
      if (obs->equiv(vbs))
        return true; // OBS=VBS RHF, UHF
      else
        return false; // VBS != OBS
    }
  }
  return false; // general reference
}

const Ref<OrbitalSpace>&
R12IntEval::hj_x_P(SpinCase1 S)
{
    switch(r12world()->r12tech()->ansatz()->orbital_product_GG()) {
    case R12Technology::OrbProdGG_ij:
	return(hj_i_P(S));
    case R12Technology::OrbProdGG_pq:
	return(hj_p_P(S));
    default:
	throw ProgrammingError("R12IntEval::GGspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}

const Ref<OrbitalSpace>&
R12IntEval::hj_x_p(SpinCase1 S)
{
    switch(r12world()->r12tech()->ansatz()->orbital_product_GG()) {
    case R12Technology::OrbProdGG_ij:
	return(hj_i_p(S));
    case R12Technology::OrbProdGG_pq:
	return(hj_p_p(S));
    default:
	throw ProgrammingError("R12IntEval::GGspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}

const Ref<OrbitalSpace>&
R12IntEval::hj_x_m(SpinCase1 S)
{
    switch(r12world()->r12tech()->ansatz()->orbital_product_GG()) {
    case R12Technology::OrbProdGG_ij:
        return(hj_i_m(S));
    case R12Technology::OrbProdGG_pq:
        return(hj_p_m(S));
    default:
        throw ProgrammingError("R12IntEval::GGspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}

const Ref<OrbitalSpace>&
R12IntEval::hj_x_a(SpinCase1 S)
{
    switch(r12world()->r12tech()->ansatz()->orbital_product_GG()) {
    case R12Technology::OrbProdGG_ij:
        return(hj_i_a(S));
    case R12Technology::OrbProdGG_pq:
        return(hj_p_a(S));
    default:
        throw ProgrammingError("R12IntEval::GGspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}

const Ref<OrbitalSpace>&
R12IntEval::hj_x_A(SpinCase1 S)
{
    switch(r12world()->r12tech()->ansatz()->orbital_product_GG()) {
    case R12Technology::OrbProdGG_ij:
	return(hj_i_A(S));
    case R12Technology::OrbProdGG_pq:
	return(hj_p_A(S));
    default:
	throw ProgrammingError("R12IntEval::GGspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}

const Ref<OrbitalSpace>&
R12IntEval::K_x_P(SpinCase1 S)
{
    switch(r12world()->r12tech()->ansatz()->orbital_product_GG()) {
    case R12Technology::OrbProdGG_ij:
	return(K_i_P(S));
    case R12Technology::OrbProdGG_pq:
	return(K_p_P(S));
    default:
	throw ProgrammingError("R12IntEval::GGspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}

const Ref<OrbitalSpace>&
R12IntEval::K_x_p(SpinCase1 S)
{
    switch(r12world()->r12tech()->ansatz()->orbital_product_GG()) {
    case R12Technology::OrbProdGG_ij:
	return(K_i_p(S));
    case R12Technology::OrbProdGG_pq:
	return(K_p_p(S));
    default:
	throw ProgrammingError("R12IntEval::GGspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}

const Ref<OrbitalSpace>&
R12IntEval::K_x_m(SpinCase1 S)
{
    switch(r12world()->r12tech()->ansatz()->orbital_product_GG()) {
    case R12Technology::OrbProdGG_ij:
	return(K_i_m(S));
    case R12Technology::OrbProdGG_pq:
	return(K_p_m(S));
    default:
	throw ProgrammingError("R12IntEval::GGspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}

const Ref<OrbitalSpace>&
R12IntEval::K_x_a(SpinCase1 S)
{
    switch(r12world()->r12tech()->ansatz()->orbital_product_GG()) {
    case R12Technology::OrbProdGG_ij:
	return(K_i_a(S));
    case R12Technology::OrbProdGG_pq:
	return(K_p_a(S));
    default:
	throw ProgrammingError("R12IntEval::GGspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}

const Ref<OrbitalSpace>&
R12IntEval::K_x_A(SpinCase1 S)
{
    switch(r12world()->r12tech()->ansatz()->orbital_product_GG()) {
    case R12Technology::OrbProdGG_ij:
	return(K_i_A(S));
    case R12Technology::OrbProdGG_pq:
	return(K_p_A(S));
    default:
	throw ProgrammingError("R12IntEval::GGspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}

const Ref<OrbitalSpace>&
R12IntEval::F_x_A(SpinCase1 S)
{
    switch(r12world()->r12tech()->ansatz()->orbital_product_GG()) {
    case R12Technology::OrbProdGG_ij:
	return(F_i_A(S));
    case R12Technology::OrbProdGG_pq:
	return(F_p_A(S));
    default:
	throw ProgrammingError("R12IntEval::GGspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}

const Ref<OrbitalSpace>&
R12IntEval::F_x_p(SpinCase1 S)
{
    switch(r12world()->r12tech()->ansatz()->orbital_product_GG()) {
    case R12Technology::OrbProdGG_ij:
    return(F_i_p(S));
    case R12Technology::OrbProdGG_pq:
    return(F_p_p(S));
    default:
    throw ProgrammingError("R12IntEval::GGspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}

const Ref<OrbitalSpace>&
R12IntEval::F_x_m(SpinCase1 S)
{
    switch(r12world()->r12tech()->ansatz()->orbital_product_GG()) {
    case R12Technology::OrbProdGG_ij:
    return(F_i_m(S));
    case R12Technology::OrbProdGG_pq:
    return(F_p_m(S));
    default:
    throw ProgrammingError("R12IntEval::GGspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}

const Ref<OrbitalSpace>&
R12IntEval::F_x_a(SpinCase1 S)
{
    switch(r12world()->r12tech()->ansatz()->orbital_product_GG()) {
    case R12Technology::OrbProdGG_ij:
    return(F_i_a(S));
    case R12Technology::OrbProdGG_pq:
    return(F_p_a(S));
    default:
    throw ProgrammingError("R12IntEval::GGspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}

const Ref<OrbitalSpace>&
R12IntEval::F_x_P(SpinCase1 S)
{
    switch(r12world()->r12tech()->ansatz()->orbital_product_GG()) {
    case R12Technology::OrbProdGG_ij:
    return(F_i_P(S));
    case R12Technology::OrbProdGG_pq:
    abort();
    default:
    throw ProgrammingError("R12IntEval::GGspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}

const Ref<OrbitalSpace>&
R12IntEval::J_x_p(SpinCase1 S)
{
    switch(r12world()->r12tech()->ansatz()->orbital_product_GG()) {
    case R12Technology::OrbProdGG_ij:
    return(J_i_p(S));
    case R12Technology::OrbProdGG_pq:
      abort();
    default:
    throw ProgrammingError("R12IntEval::GGspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}



namespace {
  /// Convert 2 spaces to SpinCase2
    SpinCase2
    spincase2(const Ref<OrbitalSpace>& space1,
              const Ref<OrbitalSpace>& space2)
    {
      char id1 = space1->id()[0];
      char id2 = space2->id()[0];
      if (id1 < 'a' && id2 < 'a')
        return AlphaAlpha;
      if (id1 < 'a' && id2 >= 'a')
        return AlphaBeta;
      if (id1 >= 'a' && id2 >= 'a')
        return BetaBeta;
      throw ProgrammingError("spincase2(space1,space2) -- BetaAlpha spaces are not allowed",
                             __FILE__,__LINE__);
    }
    std::string
    id(SpinCase2 S) {
      switch(S) {
        case AlphaBeta:  return "ab";
        case AlphaAlpha:  return "aa";
        case BetaBeta:  return "bb";
        default: MPQC_ASSERT(false); // unreachable
      }
      return std::string(); // dummy return statement to pacify picky compilers
    }
};

std::string
R12IntEval::transform_label(const Ref<OrbitalSpace>& space1,
                            const Ref<OrbitalSpace>& space2,
                            const Ref<OrbitalSpace>& space3,
                            const Ref<OrbitalSpace>& space4,
                            const std::string& operator_label) const
{
  std::ostringstream oss;
  // use physicists' notation
  oss << "<" << space1->id() << " " << space3->id() << (operator_label.empty() ? "|" : operator_label)
    << space2->id() << " " << space4->id() << ">";
  // for case-insensitive file systems append spincase
  oss << "_" << id(spincase2(space1,space3));
  return oss.str();
}

std::string
R12IntEval::transform_label(const Ref<OrbitalSpace>& space1,
                            const Ref<OrbitalSpace>& space2,
                            const Ref<OrbitalSpace>& space3,
                            const Ref<OrbitalSpace>& space4,
                            unsigned int f12,
                            const std::string& operator_label) const
{
  std::ostringstream oss;
  // use physicists' notation
  oss << "<" << space1->id() << " " << space3->id() << "| " << (operator_label.empty() ? corrfactor()->label() : operator_label)
      << "[" << f12 << "] |" << space2->id() << " " << space4->id() << ">";
  // for case-insensitive file systems append spincase
  oss << "_" << id(spincase2(space1,space3));
  return oss.str();
}

std::string
R12IntEval::transform_label(const Ref<OrbitalSpace>& space1,
                            const Ref<OrbitalSpace>& space2,
                            const Ref<OrbitalSpace>& space3,
                            const Ref<OrbitalSpace>& space4,
                            unsigned int f12_left,
                            unsigned int f12_right,
                            const std::string& operator_label) const
{
  std::ostringstream oss;
  // use physicists' notation
  oss << "<" << space1->id() << " " << space3->id() << "| " << (operator_label.empty() ? corrfactor()->label() : operator_label)
      << "[" << f12_left << "," << f12_right << "] |" << space2->id()
      << " " << space4->id() << ">";
  // for case-insensitive file systems append spincase
  oss << "_" << id(spincase2(space1,space3));
  return oss.str();
}

void
R12IntEval::spinadapt_mospace_labels(SpinCase1 spin, std::string& id, std::string& name) const
{
  // do nothing if spin-restricted
  if (!spin_polarized())
    return;

  // Prepend spin case to name
  name = prepend_spincase(spin,name);
  // Convert all characters in id which appear before '_' or '(' to upper case, if Alpha
  if (spin == Alpha) {
    std::string::const_iterator end = id.end();
    for(std::string::iterator c=id.begin(); c!=end; c++) {
      if (*c == '_' || *c == '(')
        return;
      if (*c > 'A' && *c < 'Z')
        throw ProgrammingError("R12IntEval::spinadapt() -- id should be all lower-case characters before '_'",__FILE__,__LINE__);
      if (*c > 'a' && *c < 'z') {
        *c -= 'a' - 'A';
      }
    }
  }
}



/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
