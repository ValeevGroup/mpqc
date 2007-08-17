//
// r12int_eval.cc
//
// Copyright (C) 2004 Edward Valeev
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

#ifdef __GNUG__
#pragma implementation
#endif

#include <util/misc/formio.h>
#include <util/class/scexception.h>
#include <util/ref/ref.h>
#include <util/state/state_bin.h>
#include <math/scmat/local.h>
#include <chemistry/qc/wfn/wfn.h>
#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>
#include <chemistry/qc/mbptr12/pairiter.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/transform_factory.h>
#include <chemistry/qc/mbptr12/utils.h>
#include <chemistry/qc/mbptr12/r12_amps.h>
#include <chemistry/qc/mbptr12/compute_tbint_tensor.h>
#include <chemistry/qc/mbptr12/contract_tbint_tensor.h>
#include <chemistry/qc/mbptr12/container.h>
#include <chemistry/qc/mbptr12/creator.h>

using namespace std;
using namespace sc;

#define TEST_FOCK 0
#define NOT_INCLUDE_DIAGONAL_VXB_CONTIBUTIONS 0
#define INCLUDE_EBC_CODE 1
#define INCLUDE_GBC_CODE 1
#define USE_TENSOR_CODE 1

inline int max(int a,int b) { return (a > b) ? a : b;}

// Remove when gcc has gotten rid of the bug which causes missing ERI_to_T2::I2T body when
// including definition from compute_tbint_tensor.h
#if 0
    /// Tensor elements are <pq||rs>
    class I_to_T {
    public:
      static double I2T(double I, int i1, int i3, int i2, int i4,
      const RefDiagSCMatrix& evals1,
      const RefDiagSCMatrix& evals2,
      const RefDiagSCMatrix& evals3,
      const RefDiagSCMatrix& evals4)
      {
        return I;
      }
    };
    
    /// MP2 T2 tensor elements are <ij||ab> /(e_i + e_j - e_a - e_b)
    class ERI_to_T2 {
    public:
      static double I2T(double I, int i1, int i3, int i2, int i4,
      const RefDiagSCMatrix& evals1,
      const RefDiagSCMatrix& evals2,
      const RefDiagSCMatrix& evals3,
      const RefDiagSCMatrix& evals4)
      {
        const double denom = 1.0/(evals1(i1) + evals3(i3) - evals2(i2) - evals4(i4));
        return I*denom;
      }
    };
    
    /// MP2 pseudo-T2 (S2) tensor elements are <ij||ab> /sqrt(|e_i + e_j - e_a - e_b|) such
    /// that MP2 pair energies are the diagonal elements of S2 * S2.t()
    class ERI_to_S2 {
    public:
      static double I2T(double I, int i1, int i3, int i2, int i4,
                 const RefDiagSCMatrix& evals1,
                 const RefDiagSCMatrix& evals2,
                 const RefDiagSCMatrix& evals3,
                 const RefDiagSCMatrix& evals4)
      {
        const double denom = 1.0/sqrt(fabs(evals2(i2) + evals4(i4) - evals1(i1) - evals3(i3)));
        return I*denom;
      }
    };
#endif

/*-----------------
  R12IntEval
 -----------------*/
static ClassDesc R12IntEval_cd(
  typeid(R12IntEval),"R12IntEval",2,"virtual public SavableState",
  0, 0, 0);

R12IntEval::R12IntEval(const Ref<R12IntEvalInfo>& r12i) :
  r12info_(r12i), evaluated_(false), debug_(0)
{
  this->reference();   // increase count so that I can safely create and destroy Ref<> to this
  int naocc_a, naocc_b;
  int navir_a, navir_b;
  int nall_a, nall_b;
  if (!spin_polarized()) {
    const int nocc_act = r12info_->refinfo()->docc_act()->rank();
    const int nvir_act = r12info_->vir_act()->rank();
    const int nall = r12info_->refinfo()->orbs(Alpha)->rank();
    naocc_a = naocc_b = nocc_act;
    navir_a = navir_b = nvir_act;
    nall_a = nall_b = nall;
  }
  else {
    naocc_a = occ_act(Alpha)->rank();
    naocc_b = occ_act(Beta)->rank();
    navir_a = vir_act(Alpha)->rank();
    navir_b = vir_act(Beta)->rank();
    nall_a = r12info_->refinfo()->orbs(Alpha)->rank();
    nall_b = r12info_->refinfo()->orbs(Beta)->rank();
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

  switch(r12info()->ansatz()->orbital_product()) {
  case LinearR12::OrbProd_ij:
      dim_xy_[AlphaAlpha] = dim_oo_[AlphaAlpha];
      dim_xy_[AlphaBeta] = dim_oo_[AlphaBeta];
      dim_xy_[BetaBeta] = dim_oo_[BetaBeta];
      break;
  case LinearR12::OrbProd_pq:
      {
        const unsigned int norbs_a = r12info()->refinfo()->orbs(Alpha)->rank();
        const unsigned int norbs_b = r12info()->refinfo()->orbs(Beta)->rank();
        dim_xy_[AlphaAlpha] = new SCDimension((norbs_a*(norbs_a-1))/2);
        dim_xy_[AlphaBeta] = new SCDimension(norbs_a*norbs_b);
        dim_xy_[BetaBeta] = new SCDimension((norbs_b*(norbs_b-1))/2);
      }
      break;
  default:
      throw ProgrammingError("R12IntEval::R12IntEval -- invalid orbital_product for the R12 ansatz",__FILE__,__LINE__);
  }
  for(int s=0; s<NSpinCases2; s++) {
    dim_f12_[s] = new SCDimension(corrfactor()->nfunctions()*dim_xy_[s].n());
  }
  
  if (!spin_polarized()) {
    dim_ij_s_ = new SCDimension((naocc_a*(naocc_a+1))/2);
    dim_ij_t_ = new SCDimension((naocc_a*(naocc_a-1))/2);
  }

  Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();
  for(int s=0; s<NSpinCases2; s++) {
    if (spin_polarized() || s != BetaBeta) {
      V_[s] = local_matrix_kit->matrix(dim_f12_[s],dim_oo_[s]);
      X_[s] = local_matrix_kit->matrix(dim_f12_[s],dim_f12_[s]);
      B_[s] = local_matrix_kit->matrix(dim_f12_[s],dim_f12_[s]);
      if (stdapprox() == LinearR12::StdApprox_B)
        BB_[s] = local_matrix_kit->matrix(dim_f12_[s],dim_f12_[s]);
      if (stdapprox() == LinearR12::StdApprox_C)
        BC_[s] = local_matrix_kit->matrix(dim_f12_[s],dim_f12_[s]);
      if (ebc() == false) {
        A_[s] = local_matrix_kit->matrix(dim_f12_[s],dim_vv_[s]);
        if (ks_ebcfree()) {
          Ac_[s] = local_matrix_kit->matrix(dim_f12_[s],dim_vv_[s]);
        }
#if 0
        T2_[s] = local_matrix_kit->matrix(dim_oo_[s],dim_vv_[s]);
        F12_[s] = local_matrix_kit->matrix(dim_f12_[s],dim_vv_[s]);
#endif
      }
      emp2pair_[s] = local_matrix_kit->vector(dim_oo_[s]);
    }
    else {
      V_[BetaBeta] = V_[AlphaAlpha];
      X_[BetaBeta] = X_[AlphaAlpha];
      B_[BetaBeta] = B_[AlphaAlpha];
      BB_[BetaBeta] = BB_[AlphaAlpha];
      BC_[BetaBeta] = BC_[AlphaAlpha];
      A_[BetaBeta] = A_[AlphaAlpha];
      Ac_[BetaBeta] = Ac_[AlphaAlpha];
#if 0
      T2_[BetaBeta] = T2_[AlphaAlpha];
      F12_[BetaBeta] = F12_[AlphaAlpha];
#endif
      emp2pair_[BetaBeta] = emp2pair_[AlphaAlpha];
    }
  }
  
  // WARNING Can use R12IntsAcc_MemoryGrp only for MP2 computations
  bool mp2_only = false;
  {
    Ref<LinearR12::NullCorrelationFactor> nullptr; nullptr << r12info()->corrfactor();
    if (nullptr.nonnull())
      mp2_only = true;
  }
  if (r12info()->ints_method() != R12IntEvalInfo::StoreMethod::posix && !mp2_only)
    throw InputError("R12IntEval::R12IntEval() -- the only supported storage method is posix");
  
  init_tforms_();
  // init_intermeds_ may require initialized transforms
  init_intermeds_();

  // compute canonical space of virtuals if VBS != OBS
  //if (r12info()->basis_vir()->equiv(r12info()->basis())) {
    form_canonvir_space_();
  //}
  
  Amps_ = new F12Amplitudes(this);

  this->dereference();   // to match reference() above
}

R12IntEval::R12IntEval(StateIn& si) : SavableState(si)
{
  r12info_ << SavableState::restore_state(si);

  Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();
  for(int s=0; s<NSpinCases2; s++) {
    dim_oo_[s] << SavableState::restore_state(si);
    dim_vv_[s] << SavableState::restore_state(si);
    dim_f12_[s] << SavableState::restore_state(si);
    if (si.version(::class_desc<R12IntEval>()) >= 2) {
	dim_xy_[s] << SavableState::restore_state(si);
    }
    if (!(spin_polarized() && s == BetaBeta)) {
      V_[s] = local_matrix_kit->matrix(dim_f12_[s],dim_oo_[s]);
      X_[s] = local_matrix_kit->matrix(dim_f12_[s],dim_f12_[s]);
      B_[s] = local_matrix_kit->matrix(dim_f12_[s],dim_f12_[s]);
      if (stdapprox() == LinearR12::StdApprox_B)
        BB_[s] = local_matrix_kit->matrix(dim_f12_[s],dim_f12_[s]);
      if (stdapprox() == LinearR12::StdApprox_C)
        BC_[s] = local_matrix_kit->matrix(dim_f12_[s],dim_f12_[s]);
      if (ebc() == false) {
        A_[s] = local_matrix_kit->matrix(dim_f12_[s],dim_vv_[s]);
        if (ks_ebcfree()) {
          Ac_[s] = local_matrix_kit->matrix(dim_f12_[s],dim_vv_[s]);
        }
#if 0
        T2_[s] = local_matrix_kit->matrix(dim_oo_[s],dim_vv_[s]);
        F12_[s] = local_matrix_kit->matrix(dim_f12_[s],dim_vv_[s]);
#endif
      }
      emp2pair_[s] = local_matrix_kit->vector(dim_vv_[s]);
      
      V_[s].restore(si);
      X_[s].restore(si);
      B_[s].restore(si);
      BB_[s].restore(si);
      BC_[s].restore(si);
      A_[s].restore(si);
      Ac_[s].restore(si);
#if 0
      T2_[s].restore(si);
      F12_[s].restore(si);
#endif
      emp2pair_[s].restore(si);
    }
    else {
      V_[BetaBeta] = V_[AlphaAlpha];
      X_[BetaBeta] = X_[AlphaAlpha];
      B_[BetaBeta] = B_[AlphaAlpha];
      BB_[BetaBeta] = BB_[AlphaAlpha];
      A_[BetaBeta] = A_[AlphaAlpha];
      Ac_[BetaBeta] = Ac_[AlphaAlpha];
#if 0
      T2_[BetaBeta] = T2_[AlphaAlpha];
      F12_[BetaBeta] = F12_[AlphaAlpha];
#endif
      emp2pair_[BetaBeta] = emp2pair_[AlphaAlpha];
    }
  }
  
  int num_tforms;
  si.get(num_tforms);
  for(int t=0; t<num_tforms; t++) {
    std::string tform_name;
    si.get(tform_name);
    Ref<TwoBodyMOIntsTransform> tform;
    tform << SavableState::restore_state(si);
    tform_map_[tform_name] = tform;
  }

  int evaluated; si.get(evaluated); evaluated_ = (bool) evaluated;
  si.get(debug_);

  init_tforms_();
}

R12IntEval::~R12IntEval()
{
}

void
R12IntEval::save_data_state(StateOut& so)
{
  SavableState::save_state(r12info_.pointer(),so);

  for(int s=0; s<NSpinCases2; s++) {
    SavableState::save_state(dim_oo_[s].pointer(),so);
    SavableState::save_state(dim_vv_[s].pointer(),so);
    SavableState::save_state(dim_f12_[s].pointer(),so);
    SavableState::save_state(dim_xy_[s].pointer(),so);
    if (!(spin_polarized() && s == BetaBeta)) {
      V_[s].save(so);
      X_[s].save(so);
      B_[s].save(so);
      BB_[s].save(so);
      BC_[s].save(so);
      A_[s].save(so);
      Ac_[s].save(so);
#if 0
      T2_[s].save(so);
      F12_[s].save(so);
#endif
      emp2pair_[s].save(so);
    }
  }

  int num_tforms = tform_map_.size();
  so.put(num_tforms);
  TformMap::iterator first_tform = tform_map_.begin();
  TformMap::iterator last_tform = tform_map_.end();
  for(TformMap::iterator t=first_tform; t!=last_tform; t++) {
    so.put((*t).first);
    SavableState::save_state((*t).second.pointer(),so);
  }

  so.put((int)evaluated_);
  so.put(debug_);
}

void
R12IntEval::obsolete()
{
  evaluated_ = false;

  // make all transforms obsolete
  TformMap::iterator first_tform = tform_map_.begin();
  TformMap::iterator last_tform = tform_map_.end();
  for(TformMap::iterator t=first_tform; t!=last_tform; t++) {
    (*t).second->obsolete();
  }

  init_intermeds_();
}

void R12IntEval::set_debug(int debug) { if (debug >= 0) { debug_ = debug; r12info_->set_debug_level(debug_); }};
void R12IntEval::set_dynamic(bool dynamic) { r12info_->set_dynamic(dynamic); };
void R12IntEval::set_print_percent(double pp) { r12info_->set_print_percent(pp); };
void R12IntEval::set_memory(size_t nbytes) { r12info_->set_memory(nbytes); };

const Ref<R12IntEvalInfo>& R12IntEval::r12info() const { return r12info_; };
RefSCDimension R12IntEval::dim_oo_s() const { return dim_ij_s_; };
RefSCDimension R12IntEval::dim_oo_t() const { return dim_ij_t_; };
RefSCDimension R12IntEval::dim_oo(SpinCase2 S) const { return dim_oo_[S]; }
RefSCDimension R12IntEval::dim_vv(SpinCase2 S) const { return dim_vv_[S]; }
RefSCDimension R12IntEval::dim_f12(SpinCase2 S) const { return dim_f12_[S]; }
RefSCDimension R12IntEval::dim_xy(SpinCase2 S) const { return dim_xy_[S]; }

const RefSCMatrix&
R12IntEval::V(SpinCase2 S) {
  compute();
  if (!spin_polarized() && (S == AlphaAlpha || S == BetaBeta))
    antisymmetrize(V_[AlphaAlpha],V_[AlphaBeta],
                   xspace(Alpha),
                   occ_act(Alpha));
  return V_[S];
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
                   xspace(Alpha),
                   xspace(Alpha));
  return to_lower_triangle(X_[S]);
}

RefSymmSCMatrix
R12IntEval::B(SpinCase2 S) {
  compute();
  if (!spin_polarized() && (S == AlphaAlpha || S == BetaBeta))
    antisymmetrize(B_[AlphaAlpha],B_[AlphaBeta],
                   xspace(Alpha),
                   xspace(Alpha));
  return to_lower_triangle(B_[S]);
}

RefSymmSCMatrix
R12IntEval::BB(SpinCase2 S) {
  if (stdapprox() != LinearR12::StdApprox_B)
    throw ProgrammingError("R12IntEval::BB() -- called but standard approximation is not B",__FILE__,__LINE__);
  compute();
  if (!spin_polarized() && (S == AlphaAlpha || S == BetaBeta))
    antisymmetrize(BB_[AlphaAlpha],BB_[AlphaBeta],
                   xspace(Alpha),
                   xspace(Alpha));
  return to_lower_triangle(BB_[S]);
}

RefSymmSCMatrix
R12IntEval::BC(SpinCase2 S) {
  if (stdapprox() != LinearR12::StdApprox_C)
    throw ProgrammingError("R12IntEval::BC() -- called but standard approximation is not C",__FILE__,__LINE__);
  compute();
  if (!spin_polarized() && (S == AlphaAlpha || S == BetaBeta))
    antisymmetrize(BC_[AlphaAlpha],BC_[AlphaBeta],
                   xspace(Alpha),
                   xspace(Alpha));
  return to_lower_triangle(BC_[S]);
}

const RefSCMatrix&
R12IntEval::A(SpinCase2 S) {
  compute();
  if (!spin_polarized() && (S == AlphaAlpha || S == BetaBeta))
    antisymmetrize(A_[AlphaAlpha],A_[AlphaBeta],
                   xspace(Alpha),
                   vir_act(Alpha));
  return A_[S];
}

const RefSCMatrix&
R12IntEval::Ac(SpinCase2 S) {
  compute();
  if (!spin_polarized() && (S == AlphaAlpha || S == BetaBeta))
    antisymmetrize(Ac_[AlphaAlpha],Ac_[AlphaBeta],
                   xspace(Alpha),
                   vir_act(Alpha));
  return Ac_[S];
}

const RefSCMatrix&
R12IntEval::T2(SpinCase2 S) {
  compute();
#if 0
  if (!spin_polarized() && (S == AlphaAlpha || S == BetaBeta))
    antisymmetrize(T2_[AlphaAlpha],T2_[AlphaBeta],
                   occ_act(Alpha),
                   vir_act(Alpha));
  return T2_[S];
#else
  return amps()->T2(S);
#endif
}

const RefSCMatrix&
R12IntEval::F12(SpinCase2 S) {
  compute();
#if 0
  if (!spin_polarized() && (S == AlphaAlpha || S == BetaBeta))
    antisymmetrize(F12_[AlphaAlpha],F12_[AlphaBeta],
                   xspace(Alpha),
                   vir_act(Alpha));
  return F12_[S];
#else
  return amps()->Fvv(S);
#endif
}

Ref<F12Amplitudes>
R12IntEval::amps()
{
  return Amps_;
}

const RefSCVector&
R12IntEval::emp2(SpinCase2 S)
{
  compute();
  if (!spin_polarized() && S == BetaBeta)
    return emp2pair_[AlphaAlpha];
  else
    return emp2pair_[S];
}

const RefDiagSCMatrix&
R12IntEval::evals(SpinCase1 S) const {
  if (spin_polarized()) {
    if (S == Alpha)
      return r12info()->refinfo()->orbs(Alpha)->evals();
    else
      return r12info()->refinfo()->orbs(Beta)->evals();
  }
  else
    return r12info_->refinfo()->orbs()->evals();
}

RefDiagSCMatrix R12IntEval::evals() const {
  if (spin_polarized())
    throw ProgrammingError("R12IntEval::evals() called but reference determinant spin-polarized",
                           __FILE__,__LINE__);

  return r12info_->refinfo()->orbs()->evals();
};

RefDiagSCMatrix R12IntEval::evals_a() const {
  return r12info()->refinfo()->orbs(Alpha)->evals();
}

RefDiagSCMatrix R12IntEval::evals_b() const {
  return r12info()->refinfo()->orbs(Beta)->evals();
}

void
R12IntEval::checkpoint_() const
{
  int me = r12info_->msg()->me();
  Wavefunction* wfn = r12info_->wfn();

  if (me == 0 && wfn->if_to_checkpoint()) {
    StateOutBin stateout(wfn->checkpoint_file());
    SavableState::save_state(wfn,stateout);
    ExEnv::out0() << indent << "Checkpointed Wavefunction" << endl;
  }
}

void
R12IntEval::init_tforms_()
{
  // Should be moved to Transform manager
}

Ref<TwoBodyMOIntsTransform>
R12IntEval::get_tform_(const std::string& tform_name) const
{
  TformMap::const_iterator tform_iter = tform_map_.find(tform_name);
  if (tform_iter == tform_map_.end()) {
    std::string errmsg = "R12IntEval::get_tform_() -- transform " + tform_name + " is not known";
    throw TransformNotFound(errmsg.c_str(),__FILE__,__LINE__);
  }

  return (*tform_iter).second;
}

void
R12IntEval::add_tform(const std::string& label,
                      const Ref<TwoBodyMOIntsTransform>& T)
{
  tform_map_[label] = T;
}

void
R12IntEval::init_intermeds_()
{
  for(int s=0; s<NSpinCases2; s++) {
    V_[s].assign(0.0);
    X_[s].assign(0.0);
    B_[s].assign(0.0);
    if (stdapprox() == LinearR12::StdApprox_B)
      BB_[s].assign(0.0);
    if (stdapprox() == LinearR12::StdApprox_C)
      BC_[s].assign(0.0);
    emp2pair_[s].assign(0.0);
    if (ebc() == false) {
      A_[s].assign(0.0);
      if (ks_ebcfree()) {
        Ac_[s].assign(0.0);
      }
#if 0
      T2_[s].assign(0.0);
      F12_[s].assign(0.0);
#endif
    }
  }
  
  // nothing to do if no explicit correlation
  Ref<LinearR12::NullCorrelationFactor> no12ptr; no12ptr << corrfactor();
  if (no12ptr.nonnull())
    return;
  
  Ref<LinearR12::G12CorrelationFactor> g12ptr; g12ptr << corrfactor();
  Ref<LinearR12::G12NCCorrelationFactor> g12ncptr; g12ncptr << corrfactor();
  Ref<LinearR12::GenG12CorrelationFactor> gg12ptr; gg12ptr << corrfactor();
  Ref<LinearR12::R12CorrelationFactor> r12ptr; r12ptr << corrfactor();
  if (r12ptr.nonnull()) {
    init_intermeds_r12_();
  }
  else if (g12ptr.nonnull() || g12ncptr.nonnull() || gg12ptr.nonnull()) {
    init_intermeds_g12_();
  }
}

void
R12IntEval::init_intermeds_r12_()
{
  if (r12info_->msg()->me() == 0) {

    for(int s=0; s<nspincases2(); s++) {

	const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
	const Ref<MOIndexSpace>& occ1_act = occ_act(case1(spincase2));
	const Ref<MOIndexSpace>& occ2_act = occ_act(case2(spincase2));
	const Ref<MOIndexSpace>& xspace1 = xspace(case1(spincase2));
	const Ref<MOIndexSpace>& xspace2 = xspace(case2(spincase2));

#if NOT_INCLUDE_DIAGONAL_VXB_CONTIBUTIONS
	V_[s]->assign(0.0);
	B_[s]->assign(0.0);
#else

	// compute identity operator in xc.pair/act.occ.pair basis
	RefSCMatrix I = compute_I_(xspace1,xspace2,occ1_act,occ2_act);
	if (spincase2 == AlphaBeta) {
	    V_[s].accumulate(I);
	}
	else {
	    antisymmetrize(V_[s],I,xspace1,occ1_act);
	}

	B_[s]->unit();
	if (stdapprox() == LinearR12::StdApprox_C)
	    BC_[s]->unit();
#endif
    }
  }
#if !NOT_INCLUDE_DIAGONAL_VXB_CONTIBUTIONS
  r2_contrib_to_X_new_();
#endif
}

/// Compute <space1 space2|space3 space4>
RefSCMatrix
R12IntEval::compute_I_(const Ref<MOIndexSpace>& space1,
		       const Ref<MOIndexSpace>& space2,
		       const Ref<MOIndexSpace>& space3,
		       const Ref<MOIndexSpace>& space4)
{
  /*--------------------------
    Compute overlap integrals
   --------------------------*/
  RefSCMatrix S_13;
  r12info_->compute_overlap_ints(space1, space3, S_13);
  RefSCMatrix S_24;
  if (space1 == space2 && space3 == space4) {
    S_24 = S_13;
  }
  else {
    r12info_->compute_overlap_ints(space2, space4, S_24);
  }
  const int nproc = r12info_->msg()->n();
  const int me = r12info_->msg()->me();

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

  r12info_->msg()->sum(I_array,n1234);

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
R12IntEval::compute_r2_(const Ref<MOIndexSpace>& space1,
                        const Ref<MOIndexSpace>& space2,
                        const Ref<MOIndexSpace>& space3,
                        const Ref<MOIndexSpace>& space4)
{
  /*-----------------------------------------------------
    Compute overlap, dipole, quadrupole moment integrals
   -----------------------------------------------------*/
  RefSCMatrix S_13, MX_13, MY_13, MZ_13, MXX_13, MYY_13, MZZ_13;
  r12info_->compute_multipole_ints(space1, space3, MX_13, MY_13, MZ_13, MXX_13, MYY_13, MZZ_13);
  r12info_->compute_overlap_ints(space1, space3, S_13);

  RefSCMatrix S_24, MX_24, MY_24, MZ_24, MXX_24, MYY_24, MZZ_24;
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
    r12info_->compute_multipole_ints(space2, space4, MX_24, MY_24, MZ_24, MXX_24, MYY_24, MZZ_24);
    r12info_->compute_overlap_ints(space2, space4, S_24);
  }
  if (debug_ >= DefaultPrintThresholds::diagnostics)
    ExEnv::out0() << indent << "Computed overlap and multipole moment integrals" << endl;

  const int nproc = r12info_->msg()->n();
  const int me = r12info_->msg()->me();

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

  r12info_->msg()->sum(r2_array,n1234);

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
  unsigned int me = r12info_->msg()->me();

  for(int s=0; s<nspincases2(); s++) {

    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
    const Ref<MOIndexSpace>& space1 = xspace(case1(spincase2));
    const Ref<MOIndexSpace>& space2 = xspace(case2(spincase2));
    
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


const Ref<MOIndexSpace>&
R12IntEval::focc(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return focc(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  // compute the Fock matrix between the complement and all occupieds and
  // create the new Fock-weighted space
  if (focc_space_[s].null()) {
    const Ref<MOIndexSpace>& occ_space = occ(spin);
    const Ref<MOIndexSpace>& cabs_space = r12info()->ribs_space(spin);
    
    RefSCMatrix F_ri_o = fock_(cabs_space,occ_space,spin);
    if (debug_ >= DefaultPrintThresholds::allN2)
      F_ri_o.print("Fock matrix (CABS/occ.)");

    std::string id = "m_F(a')";
    std::string name = "Fock-weighted occupied MOs sorted by energy";
    spinadapt_mospace_labels(spin,id,name);
    
    focc_space_[s] = new MOIndexSpace(id, name, occ_space, cabs_space->coefs()*F_ri_o,
                                      cabs_space->basis());
  }
  
  return focc_space_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::focc_act(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return focc_act(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  form_focc_act(spin);
  return factocc_space_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::kocc_act(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return kocc_act(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  form_focc_act(spin);
  return kactocc_space_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::hjocc_act(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return hjocc_act(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  form_focc_act(spin);
  return hjactocc_space_[s];
}

void
R12IntEval::form_focc_act(SpinCase1 spin)
{
  const unsigned int s = static_cast<unsigned int>(spin);
  // compute the Fock matrix between the complement and active occupieds and
  // create the new Fock-weighted space
  if (factocc_space_[s].null()) {
    const Ref<MOIndexSpace>& occ_space = occ(spin);
    const Ref<MOIndexSpace>& act_occ_space = r12info()->refinfo()->occ_act(spin);
    const Ref<MOIndexSpace>& cabs_space = r12info()->ribs_space(spin);
    
    RefSCMatrix FmK_ri_ao = fock_(cabs_space,act_occ_space,spin,1.0,0.0);
    RefSCMatrix K_ri_ao = exchange_(occ_space,cabs_space,act_occ_space);
    K_ri_ao.scale(-1.0);
    RefSCMatrix F_ri_ao = FmK_ri_ao.clone(); F_ri_ao.assign(FmK_ri_ao); F_ri_ao.accumulate(K_ri_ao);
    K_ri_ao.scale(-1.0);
    if (debug_ >= DefaultPrintThresholds::allN2) {
      F_ri_ao.print("Fock matrix (CABS/act.occ.)");
      K_ri_ao.print("Exchange matrix (CABS/act.occ.)");
      FmK_ri_ao.print("(h+J) matrix (CABS/act.occ.)");
    }
    
    // Fock
    {
      std::string id = "i_F(a')";
      std::string name = "Fock-weighted active occupied MOs sorted by energy";
      spinadapt_mospace_labels(spin,id,name);
      
      factocc_space_[s] = new MOIndexSpace(id, name, act_occ_space, cabs_space->coefs()*F_ri_ao,
                                           cabs_space->basis());
    }
    // Exchange
    {
      std::string id = "i_K(a')";
      std::string name = "Exchange-weighted active occupied MOs sorted by energy";
      spinadapt_mospace_labels(spin,id,name);
      
      kactocc_space_[s] = new MOIndexSpace(id, name, act_occ_space, cabs_space->coefs()*K_ri_ao,
                                           cabs_space->basis());
    }
    // Fock
    {
      std::string id = "i_hJ(a')";
      std::string name = "(h+J)-weighted active occupied MOs sorted by energy";
      spinadapt_mospace_labels(spin,id,name);
      
      hjactocc_space_[s] = new MOIndexSpace(id, name, act_occ_space, cabs_space->coefs()*FmK_ri_ao,
                                            cabs_space->basis());
    }
  }
}

const Ref<MOIndexSpace>&
R12IntEval::fvir_act(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return fvir_act(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  form_fvir_act(spin);
  return factvir_space_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::kvir_act(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return kvir_act(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  form_fvir_act(spin);
  return kactvir_space_[s];
}

void
R12IntEval::form_fvir_act(SpinCase1 spin)
{
  const unsigned int s = static_cast<unsigned int>(spin);
  // compute the Fock matrix between the complement and active virtuals and
  // create the new Fock-weighted space
  if (factvir_space_[s].null()) {
    const Ref<MOIndexSpace>& occ_space = occ(spin);
    const Ref<MOIndexSpace>& act_vir_space = vir_act(spin);
    const Ref<MOIndexSpace>& cabs_space = r12info()->ribs_space(spin);
    
    RefSCMatrix FmK_ri_av = fock_(cabs_space,act_vir_space,spin,1.0,0.0);
    RefSCMatrix K_ri_av = exchange_(occ_space,cabs_space,act_vir_space);
    K_ri_av.scale(-1.0);
    RefSCMatrix F_ri_av = FmK_ri_av; F_ri_av.accumulate(K_ri_av);
    K_ri_av.scale(-1.0);
    if (debug_ >= DefaultPrintThresholds::allN2) {
      F_ri_av.print("Fock matrix (CABS/act.vir.)");
      K_ri_av.print("Exchange matrix (CABS/act.vir.)");
    }
    
    // Fock
    {
      std::string id = "a_F(a')";
      std::string name = "Fock-weighted active virtual MOs sorted by energy";
      spinadapt_mospace_labels(spin,id,name);
      
      factvir_space_[s] = new MOIndexSpace(id, name, act_vir_space, cabs_space->coefs()*F_ri_av,
                                          cabs_space->basis());
    }
    // Exchange
    {
      std::string id = "a_K(a')";
      std::string name = "Exchange-weighted active virtual MOs sorted by energy";
      spinadapt_mospace_labels(spin,id,name);
      
      kactvir_space_[s] = new MOIndexSpace(id, name, act_vir_space, cabs_space->coefs()*K_ri_av,
                                          cabs_space->basis());
    }
  }
}

const Ref<MOIndexSpace>&
R12IntEval::kribs(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return kribs(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  form_fribs(spin);
  return kribs_space_[s];
}

void
R12IntEval::form_fribs(SpinCase1 spin)
{
  const unsigned int s = static_cast<unsigned int>(spin);
  // compute the Fock matrix between the complement and active occupieds and
  // create the new Fock-weighted space
  if (kribs_space_[s].null()) {
    const Ref<MOIndexSpace>& occ_space = occ(spin);
    const Ref<MOIndexSpace>& abs_space = r12info()->abs_space();
    const Ref<MOIndexSpace>& cabs_space = r12info()->ribs_space(spin);
    
    RefSCMatrix K_abs_ri = exchange_(occ_space,abs_space,cabs_space);
    if (debug_ >= DefaultPrintThresholds::allN2) {
      K_abs_ri.print("Exchange matrix (ABS/CABS)");
    }
    
    // Exchange
    {
      std::string id = "a'_K(p')";
      std::string name = "Exchange-weighted CABS";
      spinadapt_mospace_labels(spin,id,name);
      
      kribs_space_[s] = new MOIndexSpace(id, name, cabs_space, abs_space->coefs()*K_abs_ri,
                                         abs_space->basis());
    }
  }
}

#if 1
const Ref<MOIndexSpace>&
R12IntEval::kribs_T(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return kribs_T(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  form_fribs_T(spin);
  return kribs_T_space_[s];
}

void
R12IntEval::form_fribs_T(SpinCase1 spin)
{
  const unsigned int s = static_cast<unsigned int>(spin);
  // compute the Fock matrix between the complement and active occupieds and
  // create the new Fock-weighted space
  if (kribs_T_space_[s].null()) {
    const Ref<MOIndexSpace>& occ_space = occ(spin);
    const Ref<MOIndexSpace>& abs_space = r12info()->abs_space();
    const Ref<MOIndexSpace>& cabs_space = r12info()->ribs_space(spin);

    RefSCMatrix K_ri_abs = exchange_(occ_space,abs_space,cabs_space).t();
    if (debug_ >= DefaultPrintThresholds::allN2) {
      K_ri_abs.print("Exchange matrix (CABS/ABS)");
    }

    // Exchange
    {
      std::string id = "p'_K(a')";
      std::string name = "Exchange-weighted RIBS";
      spinadapt_mospace_labels(spin,id,name);
      
      kribs_T_space_[s] = new MOIndexSpace(id, name, abs_space, cabs_space->coefs()*K_ri_abs,
                                         cabs_space->basis());
    }
  }
}
#endif

const Ref<MOIndexSpace>&
R12IntEval::hjocc_act_obs(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return hjocc_act_obs(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  form_focc_act_obs(spin);
  return hjactocc_obs_space_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::kocc_act_obs(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return kocc_act_obs(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  form_focc_act_obs(spin);
  return kactocc_obs_space_[s];
}

void
R12IntEval::form_focc_act_obs(SpinCase1 spin)
{
  const unsigned int s = static_cast<unsigned int>(spin);
  // compute the Fock matrix between OBS and active occupieds and
  // create the new Fock-weighted space
  if (kactocc_obs_space_[s].null()) {
    const Ref<MOIndexSpace>& occ_space = occ(spin);
    const Ref<MOIndexSpace>& act_occ_space = r12info()->refinfo()->occ_act(spin);
    const Ref<MOIndexSpace>& obs_space = r12info()->refinfo()->orbs(spin);
    
    RefSCMatrix hJ_obs_ao = fock_(obs_space,act_occ_space,spin,1.0,0.0);;
    RefSCMatrix K_obs_ao = exchange_(occ_space,act_occ_space,obs_space).t();
    if (debug_ >= DefaultPrintThresholds::allN2) {
      hJ_obs_ao.print("(h+J) matrix (OBS/act.occ.)");
      K_obs_ao.print("Exchange matrix (OBS/act.occ.)");
    }

    // h+J
    {
      std::string id = "i_hJ(p)";
      std::string name = "(h+J)-weighted (through OBS) active occupied MOs sorted by energy";
      spinadapt_mospace_labels(spin,id,name);
      hjactocc_obs_space_[s] = new MOIndexSpace(id, name, act_occ_space, obs_space->coefs()*hJ_obs_ao,
						obs_space->basis());
    }

    // Exchange
    {
      std::string id = "i_K(p)";
      std::string name = "Exchange-weighted (through OBS) active occupied MOs sorted by energy";
      spinadapt_mospace_labels(spin,id,name);
      
      // as a test, use act occ space
#if 1
      kactocc_obs_space_[s] = new MOIndexSpace(id, name, act_occ_space, obs_space->coefs()*K_obs_ao,
                                           obs_space->basis());
#else
      kactocc_obs_space_[s] = new MOIndexSpace(id, name, act_occ_space, act_occ_space->coefs(),
                                              act_occ_space->basis());
#endif
    }
  }
}

const Ref<MOIndexSpace>&
R12IntEval::kvir_obs(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return kvir_obs(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  form_fvir_obs(spin);
  return kvir_obs_space_[s];
}

void
R12IntEval::form_fvir_obs(SpinCase1 spin)
{
  const unsigned int s = static_cast<unsigned int>(spin);
  // compute the Fock matrix between OBS and active occupieds and
  // create the new Fock-weighted space
  if (kvir_obs_space_[s].null()) {
    const Ref<MOIndexSpace>& occ_space = occ(spin);
    const Ref<MOIndexSpace>& vir_space = vir(spin);
    const Ref<MOIndexSpace>& obs_space = r12info()->refinfo()->orbs(spin);
    
    RefSCMatrix K_obs_vir = exchange_(occ_space,vir_space,obs_space).t();
    if (debug_ >= DefaultPrintThresholds::allN2) {
      K_obs_vir.print("Exchange matrix (OBS/vir)");
    }
    
    // Exchange
    {
      std::string id = "a_K(p)";
      std::string name = "Exchange-weighted (through OBS) virtual MOs sorted by energy";
      spinadapt_mospace_labels(spin,id,name);
      
      kvir_obs_space_[s] = new MOIndexSpace(id, name, vir_space, obs_space->coefs()*K_obs_vir,
                                           obs_space->basis());
    }
  }
}

const Ref<MOIndexSpace>&
R12IntEval::kvir_ribs(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return kvir_ribs(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  form_fvir_ribs(spin);
  return kvir_ribs_space_[s];
}

void
R12IntEval::form_fvir_ribs(SpinCase1 spin)
{
  const unsigned int s = static_cast<unsigned int>(spin);
  // compute the Fock matrix between OBS and active occupieds and
  // create the new Fock-weighted space
  if (kvir_ribs_space_[s].null()) {
    const Ref<MOIndexSpace>& occ_space = occ(spin);
    const Ref<MOIndexSpace>& vir_space = vir(spin);
    const Ref<MOIndexSpace>& ribs_space = r12info()->ribs_space();
    
    RefSCMatrix K_ribs_vir = exchange_(occ_space,vir_space,ribs_space).t();
    if (debug_ >= DefaultPrintThresholds::allN2) {
      K_ribs_vir.print("Exchange matrix (RIBS/vir)");
    }
    
    // Exchange
    {
      std::string id = "a_K(p')";
      std::string name = "Exchange-weighted (through RIBS) virtual MOs sorted by energy";
      spinadapt_mospace_labels(spin,id,name);
      
      kvir_ribs_space_[s] = new MOIndexSpace(id, name, vir_space, ribs_space->coefs()*K_ribs_vir,
                                             ribs_space->basis());
    }
  }
}

#if 1
const Ref<MOIndexSpace>&
R12IntEval::kvir_ribs_T(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return kvir_ribs_T(Alpha);

  const unsigned int s = static_cast<unsigned int>(spin);
  form_fvir_ribs_T(spin);
  return kvir_ribs_T_space_[s];
}

void
R12IntEval::form_fvir_ribs_T(SpinCase1 spin)
{
  const unsigned int s = static_cast<unsigned int>(spin);
  // compute the Fock matrix between OBS and active occupieds and
  // create the new Fock-weighted space
  if (kvir_ribs_T_space_[s].null()) {
    const Ref<MOIndexSpace>& occ_space = occ(spin);
    const Ref<MOIndexSpace>& vir_space = vir(spin);
    const Ref<MOIndexSpace>& ribs_space = r12info()->ribs_space();

    RefSCMatrix K_vir_ribs = exchange_(occ_space,vir_space,ribs_space);
    if (debug_ >= DefaultPrintThresholds::allN2) {
      K_vir_ribs.print("Exchange matrix (vir/RIBS)");
    }

      // Exchange
      {
      std::string id = "p'_K(a)";
      std::string name = "Exchange-weighted (through virtuals) RIBS";
      spinadapt_mospace_labels(spin,id,name);

      kvir_ribs_T_space_[s] = new MOIndexSpace(id, name, ribs_space, vir_space->coefs()*K_vir_ribs,
                                               vir_space->basis());
      }
    }
}

#endif

void
R12IntEval::form_canonvir_space_()
{
  // Create a complement space to all occupieds
  // Fock operator is diagonal in this space
  if (r12info_->basis_vir()->equiv(r12info_->basis())) {
    return;
  }
  
  for(int s=0; s<nspincases1(); s++) {
    const SpinCase1 spincase = static_cast<SpinCase1>(s);
    
    const Ref<MOIndexSpace>& mo_space = r12info()->refinfo()->orbs(spincase);
    const Ref<MOIndexSpace>& vir_space = r12info()->vir_sb(spincase);
    RefSCMatrix F_vir = fock_(vir_space,vir_space,spincase);
    
    int nrow = vir_space->rank();
    double* F_full = new double[nrow*nrow];
    double* F_lowtri = new double [nrow*(nrow+1)/2];
    F_vir->convert(F_full);
    int ij = 0;
    for(int row=0; row<nrow; row++) {
      int rc = row*nrow;
      for(int col=0; col<=row; col++, rc++, ij++) {
        F_lowtri[ij] = F_full[rc];
        }
      }
    RefSymmSCMatrix F_vir_lt(F_vir.rowdim(),F_vir->kit());
    F_vir_lt->assign(F_lowtri);
    F_vir = 0;
    delete[] F_full;
    delete[] F_lowtri;
    
    std::string id_sb, id;
    if (nspincases1() == 2 && spincase == Alpha) {
      id_sb = "E(sym)";
      id = "E";
    }
    else {
      id_sb = "e(sym)";
      id = "e";
    }
    Ref<MOIndexSpace> canonvir_space_symblk = new MOIndexSpace(id_sb, "VBS",
                                                               vir_space, vir_space->coefs()*F_vir_lt.eigvecs(),
                                                               vir_space->basis());
    r12info()->vir_sb(spincase, canonvir_space_symblk);
    
    RefDiagSCMatrix F_vir_evals = F_vir_lt.eigvals();
    Ref<MOIndexSpace> vir_act = new MOIndexSpace(id, "VBS",
                                                 canonvir_space_symblk->coefs(), canonvir_space_symblk->basis(),
                                                 r12info()->integral(),
                                                 F_vir_evals, 0, r12info()->refinfo()->nfzv(),
                                                 MOIndexSpace::energy);
    r12info()->vir_act(spincase, vir_act);
    Ref<MOIndexSpace> vir = new MOIndexSpace(id, "VBS",
                                             canonvir_space_symblk->coefs(), canonvir_space_symblk->basis(),
                                             r12info()->integral(),
                                             F_vir_evals, 0, 0,
                                             MOIndexSpace::energy);
    r12info()->vir(spincase, vir);
  }
}


const Ref<MOIndexSpace>&
R12IntEval::kribs_ribs(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return kribs_ribs(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  form_fribs_ribs(spin);
  return kribs_ribs_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::fribs_ribs(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return fribs_ribs(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  form_fribs_ribs(spin);
  return fribs_ribs_[s];
}

void
R12IntEval::form_fribs_ribs(SpinCase1 spin)
{
  const unsigned int s = static_cast<unsigned int>(spin);
  // compute the Fock matrix between ABS and ABS and
  // create the new Fock-weighted space
  if (fribs_ribs_[s].null()) {
    const Ref<MOIndexSpace>& occ_space = occ(spin);
    Ref<MOIndexSpace> ribs_space;
    if (r12info()->basis() != r12info()->basis_ri()) {
      ribs_space = r12info()->ribs_space();
    }
    else {
      ribs_space = r12info()->refinfo()->orbs(spin);
    }
    
    RefSCMatrix hJ_ribs_ribs = fock_(ribs_space,ribs_space,spin,1.0,0.0);
    RefSCMatrix K_ribs_ribs = exchange_(occ_space,ribs_space,ribs_space);
    K_ribs_ribs.scale(-1.0);
    RefSCMatrix F_ribs_ribs = hJ_ribs_ribs; F_ribs_ribs.accumulate(K_ribs_ribs);
    K_ribs_ribs.scale(-1.0);
    if (debug_ >= DefaultPrintThresholds::allN2) {
      hJ_ribs_ribs.print("F matrix (RIBS/RIBS)");
      K_ribs_ribs.print("Exchange matrix (RIBS/RIBS)");
    }
    
    // Fock
    {
      std::string id = "p'_F(p')";
      std::string name = "Fock-weighted (through RIBS) RIBS sorted by energy";
      spinadapt_mospace_labels(spin,id,name);
      
      fribs_ribs_[s] = new MOIndexSpace(id, name, ribs_space, ribs_space->coefs()*F_ribs_ribs,
                                            ribs_space->basis());
    }
    // Exchange
    {
      std::string id = "p'_K(p')";
      std::string name = "Exchange-weighted (through RIBS) RIBS sorted by energy";
      spinadapt_mospace_labels(spin,id,name);
      
      kribs_ribs_[s] = new MOIndexSpace(id, name, ribs_space, ribs_space->coefs()*K_ribs_ribs,
                                            ribs_space->basis());
    }
  }
}

const Ref<MOIndexSpace>&
R12IntEval::hjactocc_ribs(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return hjactocc_ribs(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  form_hjactocc_ribs(spin);
  return hjactocc_ribs_[s];
}

void
R12IntEval::form_hjactocc_ribs(SpinCase1 spin)
{
  const unsigned int s = static_cast<unsigned int>(spin);
  // compute the Fock matrix between active occupieds and RIBS and
  // create the new Fock-weighted space
  if (hjactocc_ribs_[s].null()) {
    const Ref<MOIndexSpace>& act_occ_space = occ_act(spin);
    Ref<MOIndexSpace> ribs_space;
    if (r12info()->basis() != r12info()->basis_ri()) {
      ribs_space = r12info()->ribs_space();
    }
    else {
      ribs_space = r12info()->refinfo()->orbs(spin);
    }
    
    RefSCMatrix hJ_ribs_aocc = fock_(ribs_space,act_occ_space,spin,1.0,0.0);
    if (debug_ >= DefaultPrintThresholds::allN2) {
      hJ_ribs_aocc.print("h+J matrix (RIBS/act.occ.)");
    }
    
    // Fock
    {
      std::string id = "i_hJ(p')";
      std::string name = "(h+J)-weighted (through RIBS) active occupied MOs sorted by energy";
      spinadapt_mospace_labels(spin,id,name);
      
      hjactocc_ribs_[s] = new MOIndexSpace(id, name, act_occ_space, ribs_space->coefs()*hJ_ribs_aocc,
                                           ribs_space->basis());
    }
  }
}

const Ref<MOIndexSpace>&
R12IntEval::focc_occ(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return focc_occ(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  form_focc_occ(spin);
  return focc_occ_[s];
}

void
R12IntEval::form_focc_occ(SpinCase1 spin)
{
  const unsigned int s = static_cast<unsigned int>(spin);
  // compute the Fock matrix between occupieds and occupieds and
  // create the new Fock-weighted space
  if (focc_occ_[s].null()) {
    const Ref<MOIndexSpace>& occ_space = occ(spin);
    
    RefSCMatrix F_occ = fock_(occ_space,occ_space,spin);
    if (debug_ >= DefaultPrintThresholds::allN2) {
      F_occ.print("Fock matrix (occ./occ.)");
    }
    
    // Fock
    {
      std::string id = "m_F(m)";
      std::string name = "Fock-weighted (through occupied MOs) occupied MOs sorted by energy";
      spinadapt_mospace_labels(spin,id,name);
      
      // as a test, use act occ space
      focc_occ_[s] = new MOIndexSpace(id, name, occ_space, occ_space->coefs()*F_occ,
                                      occ_space->basis());
    }
  }
}

const Ref<MOIndexSpace>&
R12IntEval::fobs_obs(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return fobs_obs(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  form_fobs_obs(spin);
  return fobs_obs_[s];
}

void
R12IntEval::form_fobs_obs(SpinCase1 spin)
{
  const unsigned int s = static_cast<unsigned int>(spin);
  // compute the Fock matrix between OBS and OBS and
  // create the new Fock-weighted space
  if (fobs_obs_[s].null()) {
    const Ref<MOIndexSpace>& obs_space = r12info()->refinfo()->orbs(spin);
    
    RefSCMatrix F_obs = fock_(obs_space,obs_space,spin);
    if (debug_ >= DefaultPrintThresholds::allN2) {
      F_obs.print("Fock matrix (OBS/OBS)");
    }
    
    // Fock
    {
      std::string id = "p_F(p)";
      std::string name = "Fock-weighted (through OBS) OBS sorted by energy";
      spinadapt_mospace_labels(spin,id,name);
      
      // as a test, use act occ space
      fobs_obs_[s] = new MOIndexSpace(id, name, obs_space, obs_space->coefs()*F_obs,
                                      obs_space->basis());
    }
  }
}

const Ref<MOIndexSpace>&
R12IntEval::focc_ribs(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return focc_ribs(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  form_focc_ribs(spin);
  return focc_ribs_[s];
}

void
R12IntEval::form_focc_ribs(SpinCase1 spin)
{
  const unsigned int s = static_cast<unsigned int>(spin);
  // compute the Fock matrix between occupieds and RIBS and
  // create the new Fock-weighted space
  if (focc_ribs_[s].null()) {
    const Ref<MOIndexSpace>& occ_space = occ(spin);
    Ref<MOIndexSpace> ribs_space;
    if (r12info()->basis() != r12info()->basis_ri()) {
      ribs_space = r12info()->ribs_space();
    }
    else {
      ribs_space = r12info()->refinfo()->orbs(spin);
    }
    
    RefSCMatrix F_ribs_occ = fock_(ribs_space,occ_space,spin);
    if (debug_ >= DefaultPrintThresholds::allN2) {
      F_ribs_occ.print("Fock matrix (RIBS/occ.)");
    }
    
    // Fock
    {
      std::string id = "m_F(p')";
      std::string name = "Fock-weighted (through RIBS) occupied MOs sorted by energy";
      spinadapt_mospace_labels(spin,id,name);
      
      focc_ribs_[s] = new MOIndexSpace(id, name, occ_space, ribs_space->coefs()*F_ribs_occ,
                                       ribs_space->basis());
    }
  }
}

const Ref<MOIndexSpace>&
R12IntEval::fobs_cabs(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return fobs_cabs(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  form_fobs_cabs(spin);
  return fobs_cabs_[s];
}

void
R12IntEval::form_fobs_cabs(SpinCase1 spin)
{
  const unsigned int s = static_cast<unsigned int>(spin);
  // compute the Fock matrix between OBS and CABS and
  // create the new Fock-weighted space
  if (fobs_cabs_[s].null()) {
    const Ref<MOIndexSpace>& obs_space = r12info()->refinfo()->orbs(spin);
    const Ref<MOIndexSpace>& cabs_space = r12info()->ribs_space(spin);
    
    RefSCMatrix F_cabs_obs = fock_(cabs_space,obs_space,spin);
    if (debug_ >= DefaultPrintThresholds::allN2) {
      F_cabs_obs.print("Fock matrix (CABS/OBS)");
    }
    
    // Fock
    {
      std::string id = "p_F(a')";
      std::string name = "Fock-weighted (through CABS) OBS sorted by energy";
      spinadapt_mospace_labels(spin,id,name);
      
      fobs_cabs_[s] = new MOIndexSpace(id, name, obs_space, cabs_space->coefs()*F_cabs_obs,
                                       cabs_space->basis());
    }
  }
}

const Ref<MOIndexSpace>&
R12IntEval::hj_i_P(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return hj_i_P(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = occ_act(spin);
  const Ref<MOIndexSpace>& intspace = r12info()->ribs_space();
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,false,true,false,
	    null,
	    hj_i_P_[s],
	    null,
	    extspace,intspace);
  return hj_i_P_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::hj_i_A(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return hj_i_A(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = occ_act(spin);
  const Ref<MOIndexSpace>& intspace = r12info()->ribs_space(spin);
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,false,true,false,
	    null,
	    hj_i_A_[s],
	    null,
	    extspace,intspace);
  return hj_i_A_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::hj_i_p(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return hj_i_p(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = occ_act(spin);
  const Ref<MOIndexSpace>& intspace = r12info()->refinfo()->orbs(spin);
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,false,true,false,
	    null,
	    hj_i_p_[s],
	    null,
	    extspace,intspace);
  return hj_i_p_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::hj_i_m(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return hj_i_m(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = occ_act(spin);
  const Ref<MOIndexSpace>& intspace = occ(spin);
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,false,true,false,
	    null,
	    hj_i_m_[s],
	    null,
	    extspace,intspace);
  return hj_i_m_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::hj_i_a(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return hj_i_a(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = occ_act(spin);
  const Ref<MOIndexSpace>& intspace = vir_act(spin);
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,false,true,false,
	    null,
	    hj_i_a_[s],
	    null,
	    extspace,intspace);
  return hj_i_a_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::hj_p_P(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return hj_p_P(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = r12info()->refinfo()->orbs(spin);
  const Ref<MOIndexSpace>& intspace = r12info()->ribs_space();
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,false,true,false,
	    null,
	    hj_p_P_[s],
	    null,
	    extspace,intspace);
  return hj_p_P_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::hj_p_A(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return hj_p_A(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = r12info()->refinfo()->orbs(spin);
  const Ref<MOIndexSpace>& intspace = r12info()->ribs_space(spin);
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,false,true,false,
	    null,
	    hj_p_A_[s],
	    null,
	    extspace,intspace);
  return hj_p_A_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::hj_p_p(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return hj_p_p(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = r12info()->refinfo()->orbs(spin);
  const Ref<MOIndexSpace>& intspace = r12info()->refinfo()->orbs(spin);
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,false,true,false,
	    null,
	    hj_p_p_[s],
	    null,
	    extspace,intspace);
  return hj_p_p_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::hj_p_m(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return hj_p_m(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = r12info()->refinfo()->orbs(spin);
  const Ref<MOIndexSpace>& intspace = occ(spin);
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,false,true,false,
	    null,
	    hj_p_m_[s],
	    null,
	    extspace,intspace);
  return hj_p_m_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::hj_p_a(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return hj_p_a(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = r12info()->refinfo()->orbs(spin);
  const Ref<MOIndexSpace>& intspace = vir_act(spin);
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,false,true,false,
	    null,
	    hj_p_a_[s],
	    null,
	    extspace,intspace);
  return hj_p_a_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::K_i_P(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_i_P(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = occ_act(spin);
  const Ref<MOIndexSpace>& intspace = r12info()->ribs_space();
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,false,false,true,
	    null,
	    null,
	    K_i_P_[s],
	    extspace,intspace);
  return K_i_P_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::K_i_A(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_i_A(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = occ_act(spin);
  const Ref<MOIndexSpace>& intspace = r12info()->ribs_space(spin);
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,false,false,true,
	    null,
	    null,
	    K_i_A_[s],
	    extspace,intspace);
  return K_i_A_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::K_i_p(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_i_p(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = occ_act(spin);
  const Ref<MOIndexSpace>& intspace = r12info()->refinfo()->orbs(spin);
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,false,false,true,
	    null,
	    null,
	    K_i_p_[s],
	    extspace,intspace);
  return K_i_p_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::K_i_m(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_i_m(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = occ_act(spin);
  const Ref<MOIndexSpace>& intspace = occ(spin);
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,false,false,true,
	    null,
	    null,
	    K_i_m_[s],
	    extspace,intspace);
  return K_i_m_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::K_i_a(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_i_a(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = occ_act(spin);
  const Ref<MOIndexSpace>& intspace = vir_act(spin);
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,false,false,true,
	    null,
	    null,
	    K_i_a_[s],
	    extspace,intspace);
  return K_i_a_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::K_m_a(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_m_a(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = occ(spin);
  const Ref<MOIndexSpace>& intspace = vir_act(spin);
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,false,false,true,
	    null,
	    null,
	    K_m_a_[s],
	    extspace,intspace);
  return K_m_a_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::K_a_a(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_a_a(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = vir_act(spin);
  const Ref<MOIndexSpace>& intspace = vir_act(spin);
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,false,false,true,
	    null,
	    null,
	    K_a_a_[s],
	    extspace,intspace);
  return K_a_a_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::K_p_P(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_p_P(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = r12info()->refinfo()->orbs(spin);
  const Ref<MOIndexSpace>& intspace = r12info()->ribs_space();
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,false,false,true,
	    null,
	    null,
	    K_p_P_[s],
	    extspace,intspace);
  return K_p_P_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::K_p_A(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_p_A(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = r12info()->refinfo()->orbs(spin);
  const Ref<MOIndexSpace>& intspace = r12info()->ribs_space(spin);
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,false,false,true,
	    null,
	    null,
	    K_p_A_[s],
	    extspace,intspace);
  return K_p_A_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::K_p_p(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_p_p(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = r12info()->refinfo()->orbs(spin);
  const Ref<MOIndexSpace>& intspace = r12info()->refinfo()->orbs(spin);
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,false,false,true,
	    null,
	    null,
	    K_p_p_[s],
	    extspace,intspace);
  return K_p_p_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::K_p_m(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_p_m(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = r12info()->refinfo()->orbs(spin);
  const Ref<MOIndexSpace>& intspace = occ(spin);
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,false,false,true,
	    null,
	    null,
	    K_p_m_[s],
	    extspace,intspace);
  return K_p_m_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::K_p_a(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_p_a(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = r12info()->refinfo()->orbs(spin);
  const Ref<MOIndexSpace>& intspace = vir_act(spin);
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,false,false,true,
	    null,
	    null,
	    K_p_a_[s],
	    extspace,intspace);
  return K_p_a_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::K_P_P(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return K_P_P(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = r12info()->ribs_space();
  const Ref<MOIndexSpace>& intspace = r12info()->ribs_space();
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,true,false,true,
	    F_P_P_[s],
	    null,
	    K_P_P_[s],
	    extspace,intspace);
  return K_P_P_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::F_P_P(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_P_P(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = r12info()->ribs_space();
  const Ref<MOIndexSpace>& intspace = r12info()->ribs_space();
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,true,false,true,
	    F_P_P_[s],
	    null,
	    K_P_P_[s],
	    extspace,intspace);
  return F_P_P_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::F_p_A(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_p_A(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = r12info()->refinfo()->orbs(spin);
  const Ref<MOIndexSpace>& intspace = r12info()->ribs_space(spin);
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,true,false,false,
	    F_p_A_[s],
	    null,
	    null,
	    extspace,intspace);
  return F_p_A_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::F_p_p(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_p_p(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = r12info()->refinfo()->orbs(spin);
  const Ref<MOIndexSpace>& intspace = r12info()->refinfo()->orbs(spin);
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,true,false,false,
	    F_p_p_[s],
	    null,
	    null,
	    extspace,intspace);
  return F_p_p_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::F_m_m(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_m_m(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = r12info()->refinfo()->occ(spin);
  const Ref<MOIndexSpace>& intspace = r12info()->refinfo()->occ(spin);
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,true,false,false,
	    F_m_m_[s],
	    null,
	    null,
	    extspace,intspace);
  return F_m_m_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::F_m_a(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_m_a(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = r12info()->refinfo()->occ(spin);
  const Ref<MOIndexSpace>& intspace = vir_act(spin);
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,true,false,false,
	    F_m_a_[s],
	    null,
	    null,
	    extspace,intspace);
  return F_m_a_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::F_m_P(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_m_P(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = r12info()->refinfo()->occ(spin);
  const Ref<MOIndexSpace>& intspace = r12info()->ribs_space();
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,true,false,false,
	    F_m_P_[s],
	    null,
	    null,
	    extspace,intspace);
  return F_m_P_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::F_m_A(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_m_A(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = r12info()->refinfo()->occ(spin);
  const Ref<MOIndexSpace>& intspace = r12info()->ribs_space(spin);
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,true,false,false,
	    F_m_A_[s],
	    null,
	    null,
	    extspace,intspace);
  return F_m_A_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::F_i_A(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_i_A(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = r12info()->refinfo()->occ_act(spin);
  const Ref<MOIndexSpace>& intspace = r12info()->ribs_space(spin);
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,true,false,false,
	    F_i_A_[s],
	    null,
	    null,
	    extspace,intspace);
  return F_i_A_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::F_a_a(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_a_a(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = vir_act(spin);
  const Ref<MOIndexSpace>& intspace = vir_act(spin);
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,true,false,false,
	    F_a_a_[s],
	    null,
	    null,
	    extspace,intspace);
  return F_a_a_[s];
}

const Ref<MOIndexSpace>&
R12IntEval::F_a_A(SpinCase1 spin)
{
  if (!spin_polarized() && spin == Beta)
    return F_a_A(Alpha);
  
  const unsigned int s = static_cast<unsigned int>(spin);
  const Ref<MOIndexSpace>& extspace = vir_act(spin);
  const Ref<MOIndexSpace>& intspace = r12info()->ribs_space(spin);
  Ref<MOIndexSpace> null;
  f_bra_ket(spin,true,false,false,
	    F_a_A_[s],
	    null,
	    null,
	    extspace,intspace);
  return F_a_A_[s];
}


void
R12IntEval::f_bra_ket(
    SpinCase1 spin,
    bool make_F,
    bool make_hJ,
    bool make_K,
    Ref<MOIndexSpace>& F,
    Ref<MOIndexSpace>& hJ,
    Ref<MOIndexSpace>& K,
    const Ref<MOIndexSpace>& extspace,
    const Ref<MOIndexSpace>& intspace
    )
{
  const unsigned int s = static_cast<unsigned int>(spin);
  const bool not_yet_computed =
      (make_F && F.null()) ||
      (make_hJ && hJ.null()) ||
      (make_K && K.null());
  if (not_yet_computed) {

      RefSCMatrix hJ_i_e;
      if (make_hJ && hJ.null()) {
	  hJ_i_e = fock_(intspace,extspace,spin,1.0,0.0);
	  if (debug_ >= DefaultPrintThresholds::allN2) {
	      std::string label("(h+J) matrix in ");
	      label += intspace->id();
	      label += "/";
	      label += extspace->id();
	      label += " basis";
	      hJ_i_e.print(label.c_str());
	  }

	  std::string id = extspace->id();  id += "_hJ(";  id += intspace->id();  id += ")";
	  std::string name = "(h+J)-weighted space";
	  name = prepend_spincase(spin,name);
	  hJ = new MOIndexSpace(id, name, extspace, intspace->coefs()*hJ_i_e,
				intspace->basis());
      }

      RefSCMatrix K_i_e;
      if (make_K && K.null()) {
	  const Ref<MOIndexSpace>& occ_space = occ(spin);
	  K_i_e = exchange_(occ_space,intspace,extspace);
	  if (debug_ >= DefaultPrintThresholds::allN2) {
	      std::string label("K matrix in ");
	      label += intspace->id();
	      label += "/";
	      label += extspace->id();
	      label += " basis";
	      K_i_e.print(label.c_str());
	  }

	  std::string id = extspace->id();  id += "_K(";  id += intspace->id();  id += ")";
	  std::string name = "K-weighted space";
	  name = prepend_spincase(spin,name);
	  K = new MOIndexSpace(id, name, extspace, intspace->coefs()*K_i_e,
				intspace->basis());
      }

      if (make_F && F.null()) {
	  RefSCMatrix F_i_e;
	  if (make_hJ) {
	      if (make_K) {
		  F_i_e = K_i_e.clone();  F_i_e.assign(K_i_e);  F_i_e.scale(-1.0);
		  F_i_e.accumulate(hJ_i_e);
	      }
	      else {
		  const Ref<MOIndexSpace>& occ_space = occ(spin);
		  F_i_e = exchange_(occ_space,intspace,extspace);
		  F_i_e.scale(-1.0);
		  F_i_e.accumulate(hJ_i_e);
	      }
	  }
	  else {
	      if (make_K) {
		  F_i_e = K_i_e.clone();  F_i_e.assign(K_i_e);  F_i_e.scale(-1.0);
		  RefSCMatrix hJ_i_e = fock_(intspace,extspace,spin,1.0,0.0);
		  F_i_e.accumulate(hJ_i_e);
	      }
	      else {
		  F_i_e = fock_(intspace,extspace,spin,1.0,1.0);
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

	  std::string id = extspace->id();  id += "_F(";  id += intspace->id();  id += ")";
	  std::string name = "F-weighted space";
          name = prepend_spincase(spin,name);
	  F = new MOIndexSpace(id, name, extspace, intspace->coefs()*F_i_e,
				intspace->basis());
      }
    
  }
}

const int
R12IntEval::tasks_with_ints_(const Ref<R12IntsAcc> ints_acc, vector<int>& map_to_twi)
{
  int nproc = r12info_->msg()->n();
  
  // Compute the number of tasks that have full access to the integrals
  // and split the work among them
  int nproc_with_ints = 0;
  for(int proc=0;proc<nproc;proc++)
    if (ints_acc->has_access(proc)) nproc_with_ints++;
 
  map_to_twi.resize(nproc);
  int count = 0;
  for(int proc=0;proc<nproc;proc++)
    if (ints_acc->has_access(proc)) {
      map_to_twi[proc] = count;
      count++;
    }
      else
        map_to_twi[proc] = -1;

  ExEnv::out0() << indent << "Computing intermediates on " << nproc_with_ints
    << " processors" << endl;
    
  return nproc_with_ints;
}

void
R12IntEval::compute()
{
  if (evaluated_)
    return;
  
  // different expressions hence codepaths depending on relationship between OBS, VBS, and RIBS
  // compare these basis sets here
  const bool obs_eq_vbs = r12info_->basis_vir()->equiv(r12info_->basis());
  const bool obs_eq_ribs = r12info_->basis_ri()->equiv(r12info_->basis());
  const LinearR12::ABSMethod absmethod = r12info()->abs_method();
  const bool cabs_method = (absmethod ==  LinearR12::ABS_CABS ||
			    absmethod == LinearR12::ABS_CABSPlus);
  // is CABS space empty? if not sure set to false and allow some crazy runtime error to happen rather than create incorrect result
  const bool cabs_empty = obs_eq_vbs && obs_eq_ribs;
  const bool vir_empty = vir(Alpha)->rank()==0 || vir(Beta)->rank()==0;
  
  Ref<LinearR12::NullCorrelationFactor> nocorrptr; nocorrptr << corrfactor();
  // if explicit correlation -- compute linear F12 theory intermediates
  if (nocorrptr.null()) {
    
    if (debug_ >= DefaultPrintThresholds::O4) {
      for(int s=0; s<nspincases2(); s++) {
        V_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"V(diag) contribution").c_str());
        X_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"X(diag) contribution").c_str());
        B_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"B(diag) contribution").c_str());
      }
    }
    
    if (obs_eq_vbs) {
#if !USE_TENSOR_CODE
      using LinearR12::TwoParticleContraction;
      using LinearR12::ABS_OBS_Contraction;
      using LinearR12::CABS_OBS_Contraction;
      for(int s=0; s<nspincases2(); s++) {
        const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
        const SpinCase1 spin1 = case1(spincase2);
        const SpinCase1 spin2 = case2(spincase2);
        
        // "ABS"-type contraction is used for projector 2 ABS/ABS+ method when OBS != RIBS
        Ref<TwoParticleContraction> tpcontract;
        if ((absmethod == LinearR12::ABS_ABS ||
             absmethod == LinearR12::ABS_ABSPlus) &&
            !obs_eq_ribs &&
            ansatz()->projector() == LinearR12::Projector_2)
          tpcontract = new ABS_OBS_Contraction(r12info()->refinfo()->orbs(spin1)->rank(),
                                               r12info()->refinfo()->occ(spin1)->rank(),
                                               r12info()->refinfo()->occ(spin2)->rank());
        else
          tpcontract = new CABS_OBS_Contraction(r12info()->refinfo()->orbs(spin1)->rank());
        contrib_to_VXB_a_new_(occ_act(spin1),
                              r12info()->refinfo()->orbs(spin1),
                              occ_act(spin2),
                              r12info()->refinfo()->orbs(spin2),
                              spincase2,tpcontract);
	if (debug_ >= DefaultPrintThresholds::O4) {
          V_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"V(diag+OBS) contribution").c_str());
          X_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"X(diag+OBS) contribution").c_str());
          B_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"B(diag+OBS) contribution").c_str());
        }
      }
      
      // These terms don't contribute to projector 3
      if (!obs_eq_ribs && ansatz()->projector() == LinearR12::Projector_2) {
        // Compute VXB using new code
        using LinearR12::Direct_Contraction;
        for(int s=0; s<nspincases2(); s++) {
          const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
          const SpinCase1 spin1 = case1(spincase2);
          const SpinCase1 spin2 = case2(spincase2);
	  Ref<MOIndexSpace> xspace1, xspace2;
	  if (cabs_method) {
	    xspace1 = r12info()->ribs_space(spin1);
	    xspace2 = r12info()->ribs_space(spin2);
	  }
	  else {
	    xspace1 = r12info()->ribs_space();
	    xspace2 = r12info()->ribs_space();
	  }
          Ref<TwoParticleContraction> tpcontract = new Direct_Contraction(
                                                         r12info()->refinfo()->occ(spin1)->rank(),
                                                         xspace2->rank(),-1.0
                                                       );
          contrib_to_VXB_a_new_(occ_act(spin1),
                                r12info()->refinfo()->occ(spin1),
                                occ_act(spin2),
                                xspace2,
                                spincase2,tpcontract);
          
          if (spincase2 == AlphaBeta && occ_act(spin1) != occ_act(spin2)) {
            Ref<TwoParticleContraction> tpcontract = new Direct_Contraction(
                                                           xspace1->rank(),
                                                           r12info()->refinfo()->occ(spin2)->rank(),-1.0
                                                         );
            contrib_to_VXB_a_new_(occ_act(spin1),
                                  xspace1,
                                  occ_act(spin2),
                                  r12info()->refinfo()->occ(spin2),
                                  spincase2,tpcontract);
          }
          
	  if (debug_ >= DefaultPrintThresholds::O4) {
            V_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"V(diag+OBS+ABS) contribution").c_str());
            X_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"X(diag+OBS+ABS) contribution").c_str());
            B_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"B(diag+OBS+ABS) contribution").c_str());
          }
        }
      }
      
#else
      contrib_to_VXB_a_();
#endif

    }
    else {
#if !USE_TENSOR_CODE
      contrib_to_VXB_gebc_vbsneqobs_();
#else
      contrib_to_VXB_a_vbsneqobs_();
#endif
    }

    // Contribution from X to B in approximation A'' is more complicated than in other methods
    // because exchange is completely skipped
    if (stdapprox() == LinearR12::StdApprox_App)
      compute_BApp_();
    
    // This is app B contribution to B -- only valid for projector 2
    if ((stdapprox() == LinearR12::StdApprox_B) &&
         ansatz()->projector() == LinearR12::Projector_2) {
      compute_BB_();
      if (debug_ >= DefaultPrintThresholds::O4)
        for(int s=0; s<nspincases2(); s++)
          BB_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"B(app. B) contribution").c_str());
    }
    
    if (stdapprox() == LinearR12::StdApprox_C) {
      compute_BC_();
      if (debug_ >= DefaultPrintThresholds::O4)
        for(int s=0; s<nspincases2(); s++)
          BC_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"B(app. C) intermediate").c_str());
    }
    
#if INCLUDE_EBC_CODE
    const bool nonzero_ebc_terms = !ebc() && !cabs_empty && !vir_empty;
    if (nonzero_ebc_terms) {
      
      // compute A, T2, and F12
      for(int s=0; s<nspincases2(); s++) {
        const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
        const SpinCase1 spin1 = case1(spincase2);
        const SpinCase1 spin2 = case2(spincase2);
        
        Ref<MOIndexSpace> xspace1 = xspace(spin1);
        Ref<MOIndexSpace> xspace2 = xspace(spin2);
        Ref<MOIndexSpace> vir1_act = vir_act(spin1);
        Ref<MOIndexSpace> vir2_act = vir_act(spin2);
        Ref<MOIndexSpace> fvir1_act = fvir_act(spin1);
        Ref<MOIndexSpace> fvir2_act = fvir_act(spin2);
        
        const Ref<SingleRefInfo> refinfo = r12info()->refinfo();
        
        std::vector<  Ref<TwoBodyMOIntsTransform> > tforms;
        Ref<R12IntEval> thisref(this);
        if (obs_eq_vbs) {
          NamedTransformCreator tform_creator(thisref,
                                              xspace1,
                                              refinfo->orbs(spin1),
                                              xspace2,
                                              refinfo->orbs(spin2),true);
          fill_container(tform_creator,tforms);
        }
        else {
          NamedTransformCreator tform_creator(thisref,
                                              xspace1,
                                              vir1_act,
                                              xspace2,
                                              vir2_act,true);
          fill_container(tform_creator,tforms);
        }
        
        compute_A_direct_(A_[s],
                          xspace1, vir1_act,
                          xspace2, vir2_act,
                          fvir1_act, fvir2_act,
                          spincase2!=AlphaBeta);
        if (ks_ebcfree()) {
          compute_A_viacomm_(Ac_[s],
                             xspace1, vir1_act,
                             xspace2, vir2_act,
                             spincase2!=AlphaBeta, tforms);
        }
      }
      
      if (debug_ >= DefaultPrintThresholds::O2N2) {
        for(int s=0; s<nspincases2(); s++) {
          const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
          std::string label = prepend_spincase(spincase2,"T2 matrix");
          amps()->T2(spincase2).print(label.c_str());
          label = prepend_spincase(spincase2,"F12(vv) matrix");
          amps()->Fvv(spincase2).print(label.c_str());
          label = prepend_spincase(spincase2,"A matrix");
          A_[s].print(label.c_str());
          if (ks_ebcfree()) {
            std::string label = prepend_spincase(spincase2,"A(comm) matrix");
            Ac_[s].print(label.c_str());
          }
        }
      }
      
      AT2_contrib_to_V_();
      // EBC contribution to B only appears in projector 2
      if (ansatz()->projector() == LinearR12::Projector_2) {
        AF12_contrib_to_B_();
      }
    }
#endif
    
#if INCLUDE_GBC_CODE
    // GBC contribution to B only appears in projector 2
    const bool nonzero_gbc_terms = !gbc() && !cabs_empty && ansatz()->projector() == LinearR12::Projector_2;
    if (nonzero_gbc_terms) {
      // These functions assume that virtuals are expanded in the same basis
      // as the occupied orbitals
      if (!obs_eq_vbs)
        throw std::runtime_error("R12IntEval::compute() -- gbc=false is only supported when basis_vir == basis");
      
      compute_B_gbc_();
    }
#endif

#if 0
    // Test new tensor compute function
    for(int s=0; s<nspincases2(); s++) {
      const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
      const SpinCase1 spin1 = case1(spincase2);
      const SpinCase1 spin2 = case2(spincase2);
      Ref<SingleRefInfo> refinfo = r12info()->refinfo();
      Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
      RefSCMatrix S2 = localkit->matrix(dim_oo(spincase2),dim_vv(spincase2));
      S2.assign(0.0);
      
      std::vector<  Ref<TwoBodyMOIntsTransform> > tforms;
      Ref<R12IntEval> thisref(this);
      NamedTransformCreator tform_creator(
        thisref,
        occ_act(spin1),
        refinfo->orbs(spin1),
        occ_act(spin2),
        refinfo->orbs(spin2)
      );
      fill_container(tform_creator,tforms);
      
      // compute S2 = <ij||ab>/sqrt(e_a+e_b-e_i-e_j) and "square" it
      compute_tbint_tensor<ManyBodyTensors::ERI_to_S2,false,false>(
        S2, corrfactor()->tbint_type_eri(),
        occ_act(spin1), vir_act(spin1),
        occ_act(spin2), vir_act(spin2),
        spincase2!=AlphaBeta, tforms
      );
      RefSCMatrix mp2pe = S2*S2.t(); mp2pe.scale(-1.0);
      mp2pe.print("S2 * S2.t : Diagonal elements should be pair energies");
      
      // another way to compute pair energies is G*T2.t()
      RefSCMatrix T2 = localkit->matrix(dim_oo(spincase2),dim_vv(spincase2));
      RefSCMatrix G = localkit->matrix(dim_oo(spincase2),dim_vv(spincase2));
      T2.assign(0.0); G.assign(0.0);
      // compute T2 and G
      compute_tbint_tensor<ManyBodyTensors::ERI_to_T2,false,false>(
        T2, corrfactor()->tbint_type_eri(),
        occ_act(spin1), vir_act(spin1),
        occ_act(spin2), vir_act(spin2),
        spincase2!=AlphaBeta, tforms
      );
      compute_tbint_tensor<ManyBodyTensors::I_to_T,false,false>(
        G, corrfactor()->tbint_type_eri(),
        occ_act(spin1), vir_act(spin1),
        occ_act(spin2), vir_act(spin2),
        spincase2!=AlphaBeta, tforms
      );
      if (debug_ >= DefaultPrintThresholds::mostO2N2)
        T2.print("T2 amplitudes");
      mp2pe = G*T2.t();
      mp2pe.print("G * T2.t : Diagonal elements should be pair energies");

      // here instead of computing G and T2 explicitly, use contract function
      {
        using namespace sc::LinearR12;
        mp2pe.assign(0.0);
        Ref<TwoParticleContraction> dircontract =
          new Direct_Contraction(vir_act(spin1)->rank(),
                                 vir_act(spin2)->rank(),
                                 1.0);
        contract_tbint_tensor<ManyBodyTensors::I_to_T,
                              ManyBodyTensors::ERI_to_T2,
                              ManyBodyTensors::I_to_T,
                              false,false,false>
          (
            mp2pe, corrfactor()->tbint_type_eri(), corrfactor()->tbint_type_eri(),
            occ_act(spin1), occ_act(spin2),
            vir_act(spin1), vir_act(spin2),
            occ_act(spin1), occ_act(spin2),
            vir_act(spin1), vir_act(spin2),
            dircontract,
            spincase2!=AlphaBeta, tforms, tforms
          );
        mp2pe.print("contract(G,T2) : Diagonal elements should be pair energies");
      }
      
      // Try computing VXB (diag+OBS) using contract function
      {
        using namespace sc::LinearR12;
        RefSCMatrix V = V_[s].clone(); V->unit();
        RefSCMatrix X = X_[s].clone(); X->assign(0.0);
        RefSCMatrix B = B_[s].clone(); B->unit();
        const Ref<SingleRefInfo> refinfo = r12info()->refinfo();
        Ref<TwoParticleContraction> tpcontract;
        const ABSMethod absmethod = r12info()->abs_method();
        if ((absmethod == LinearR12::ABS_ABS ||
          absmethod == LinearR12::ABS_ABSPlus) && !obs_eq_ribs)
          tpcontract = new ABS_OBS_Contraction(refinfo->orbs(spin1)->rank(),
                                               refinfo->occ(spin1)->rank(),
                                               refinfo->occ(spin2)->rank());
        else
          tpcontract = new CABS_OBS_Contraction(refinfo->orbs(spin1)->rank());
        
        std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_f12;
        Ref<R12IntEval> thisref(this);
        NamedTransformCreator tform_creator(
          thisref,
          occ_act(spin1),
          refinfo->orbs(spin1),
          occ_act(spin2),
          refinfo->orbs(spin2),true
        );
        fill_container(tform_creator,tforms_f12);
        
        contract_tbint_tensor<ManyBodyTensors::I_to_T,
                              ManyBodyTensors::I_to_T,
                              ManyBodyTensors::I_to_T,
                              true,false,false>
          (
            V, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_eri(),
            occ_act(spin1), occ_act(spin2),
            refinfo->orbs(spin1), refinfo->orbs(spin2),
            occ_act(spin1), occ_act(spin2),
            refinfo->orbs(spin1), refinfo->orbs(spin2),
            tpcontract,
            spincase2!=AlphaBeta, tforms_f12, tforms
          );
        V.print("contract(F12,G) + I = V(diag+OBS) if F12=R12");
        
        
        //
        // F12^2 contribution depends on the type of correlation factor
        //
        RefSCMatrix F12_sq_ijkl;
        enum {r12corrfactor, g12corrfactor, gg12corrfactor} corrfac;
        Ref<LinearR12::R12CorrelationFactor> r12corrptr; r12corrptr << corrfactor();
        Ref<LinearR12::G12CorrelationFactor> g12corrptr; g12corrptr << corrfactor();
        Ref<LinearR12::GenG12CorrelationFactor> gg12corrptr; gg12corrptr << corrfactor();
        if (r12corrptr.nonnull()) corrfac = r12corrfactor;
        if (g12corrptr.nonnull()) corrfac = g12corrfactor;
        if (gg12corrptr.nonnull()) corrfac = gg12corrfactor;
    
        switch (corrfac) {
          case r12corrfactor:  // R12^2 reduces to one-electron integrals
          F12_sq_ijkl = compute_r2_(occ_act(spin1), occ_act(spin2),
                                    occ_act(spin1), occ_act(spin2));
          break;
          
          case g12corrfactor: // G12^2 involves two-electron integrals
          case gg12corrfactor: // G12^2 involves two-electron integrals
          {
            // (i k |j l) tforms
            std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_ikjl;
            {
              NewTransformCreator tform_creator(thisref,
                                                occ_act(spin1), occ_act(spin1),
                                                occ_act(spin2), occ_act(spin2),
                                                true,true);
              fill_container(tform_creator,tforms_ikjl);
            }
            F12_sq_ijkl = X.clone();
            F12_sq_ijkl.assign(0.0);
            compute_tbint_tensor<ManyBodyTensors::I_to_T,true,true>(
              F12_sq_ijkl, corrfactor()->tbint_type_f12f12(),
              occ_act(spin1), occ_act(spin2),
              occ_act(spin1), occ_act(spin2),
              spincase2!=AlphaBeta,
              tforms_ikjl
            );
          }
          break;
          
          default:
          throw ProgrammingError("R12IntEval::compute_X_() -- unrecognized type of correlation factor",__FILE__,__LINE__);
        }
        
        if (spincase2 == AlphaAlpha || spincase2 == BetaBeta)
          antisymmetrize(X,F12_sq_ijkl,
                         occ_act(spin1),
                         occ_act(spin2));
        else
          X.accumulate(F12_sq_ijkl);

        X.print("F12^2 = X(diag)");

        contract_tbint_tensor<ManyBodyTensors::I_to_T,
                              ManyBodyTensors::I_to_T,
                              ManyBodyTensors::I_to_T,
                              true,true,false>
          (
            X, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_f12(),
            occ_act(spin1), occ_act(spin2),
            refinfo->orbs(spin1), refinfo->orbs(spin2),
            occ_act(spin1), occ_act(spin2),
            refinfo->orbs(spin1), refinfo->orbs(spin2),
            tpcontract,
            spincase2!=AlphaBeta, tforms_f12, tforms_f12
          );
        X.print("contract(F12,F12) + F12^2 = X(diag+OBS)");
        contract_tbint_tensor<ManyBodyTensors::I_to_T,
                              ManyBodyTensors::I_to_T,
                              ManyBodyTensors::I_to_T,
                              true,true,false>
          (
            B, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_t1f12(),
            occ_act(spin1), occ_act(spin2),
            refinfo->orbs(spin1), refinfo->orbs(spin2),
            occ_act(spin1), occ_act(spin2),
            refinfo->orbs(spin1), refinfo->orbs(spin2),
            tpcontract,
            spincase2!=AlphaBeta, tforms_f12, tforms_f12
          );
        contract_tbint_tensor<ManyBodyTensors::I_to_T,
                              ManyBodyTensors::I_to_T,
                              ManyBodyTensors::I_to_T,
                              true,true,false>
          (
            B, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_t2f12(),
            occ_act(spin1), occ_act(spin2),
            refinfo->orbs(spin1), refinfo->orbs(spin2),
            occ_act(spin1), occ_act(spin2),
            refinfo->orbs(spin1), refinfo->orbs(spin2),
            tpcontract,
            spincase2!=AlphaBeta, tforms_f12, tforms_f12
          );
        B.scale(0.5); RefSCMatrix Bt = B.t(); B.accumulate(Bt);
        B.print("contract(F12,[T1+T2,F12]) + I = B(diag+OBS)");
      }
      
    }
    
    // test F12 evaluator
    {
      const Ref<SingleRefInfo> refinfo = r12info()->refinfo();
      Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
      RefSCMatrix F12 = localkit->matrix(dim_f12(AlphaBeta),dim_vv(AlphaBeta));
      F12.assign(0.0);
      
      std::vector<  Ref<TwoBodyMOIntsTransform> > tforms;
      Ref<R12IntEval> thisref(this);
      NamedTransformCreator tform_creator(
        thisref,
        occ_act(Alpha),
        refinfo->orbs(Alpha),
        occ_act(Beta),
        refinfo->orbs(Beta),true
      );
      fill_container(tform_creator,tforms);
      
      compute_F12_(F12,occ_act(Alpha), vir_act(Alpha),
                   occ_act(Beta), vir_act(Beta), false, tforms);
    }
    
#endif

#if 0
    // test generic X evaluator
    {
      for(int s=0; s<nspincases2(); s++) {
        const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
        const SpinCase1 spin1 = case1(spincase2);
        const SpinCase1 spin2 = case2(spincase2);
        RefSCMatrix X;
        compute_X_(X,spincase2,occ_act(spin1),occ_act(spin2),
                               occ_act(spin1),occ_act(spin2));
        X.print("X <ii|ii>  test");
      }
    }
#endif

#if 0
    // test generic FxF evaluator
    {
      for(int s=0; s<nspincases2(); s++) {
        const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
        const SpinCase1 spin1 = case1(spincase2);
        const SpinCase1 spin2 = case2(spincase2);
        RefSCMatrix Bebc;
        compute_FxF_(Bebc,spincase2,
                     occ_act(spin1),occ_act(spin2),
                     occ_act(spin1),occ_act(spin2),
                     vir(spin1),vir(spin2),
                     vir(spin1),vir(spin2),
                     fvir_act(spin1),fvir_act(spin2));
        // -r f1 r, not r f1+f2 r computed by FxF
        Bebc.scale(-1.0);
        Bebc.print("B_{EBC} from generic FxF evaluator");
      }
    }
#endif

  }
  
  // Finally, compute MP2 energies
  const int nspincases_for_emp2pairs = (spin_polarized() ? 3 : 2);
  for(int s=0; s<nspincases_for_emp2pairs; s++) {
      const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
      const SpinCase1 spin1 = case1(spincase2);
      const SpinCase1 spin2 = case2(spincase2);
      Ref<MOIndexSpace> occ1_act = occ_act(spin1);
      Ref<MOIndexSpace> occ2_act = occ_act(spin2);
      Ref<MOIndexSpace> vir1_act = vir_act(spin1);
      Ref<MOIndexSpace> vir2_act = vir_act(spin2);
      
      std::vector<  Ref<TwoBodyMOIntsTransform> > tforms;
      Ref<R12IntEval> thisref(this);
      // If VBS==OBS and this is not a pure MP2 calculation then this tform should be available
      if (obs_eq_vbs && nocorrptr.null()) {
        NewTransformCreator tform_creator(thisref,
                                          occ1_act,
                                          r12info()->refinfo()->orbs(spin1),
                                          occ2_act,
                                          r12info()->refinfo()->orbs(spin2),1);
        fill_container(tform_creator,tforms);
      }
      else {
        NewTransformCreator tform_creator(thisref,
                                          occ1_act,
                                          vir1_act,
                                          occ2_act,
                                          vir2_act,0);
        fill_container(tform_creator,tforms);
      }
      compute_mp2_pair_energies_(emp2pair_[s],spincase2,
                                 occ1_act,vir1_act,
                                 occ2_act,vir2_act,
                                 tforms[0]);
  }
  
  // Distribute the final intermediates to every node
  globally_sum_intermeds_(true);
  
  evaluated_ = true;
  return;
}

void
R12IntEval::globally_sum_scmatrix_(RefSCMatrix& A, bool to_all_tasks, bool to_average)
{
  Ref<MessageGrp> msg = r12info_->msg();
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
  Ref<MessageGrp> msg = r12info_->msg();
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
    if (stdapprox() == LinearR12::StdApprox_B)
      globally_sum_scmatrix_(BB_[s],to_all_tasks);
    if (stdapprox() == LinearR12::StdApprox_C)
      globally_sum_scmatrix_(BC_[s],to_all_tasks);
    if (ebc() == false) {
      globally_sum_scmatrix_(A_[s],to_all_tasks);
      if (ks_ebcfree()) {
        globally_sum_scmatrix_(Ac_[s],to_all_tasks);
      }
#if 0
      globally_sum_scmatrix_(T2_[s],to_all_tasks);
      globally_sum_scmatrix_(F12_[s],to_all_tasks);
#endif
    }
    globally_sum_scvector_(emp2pair_[s],to_all_tasks);
  }

  if (debug_ >= DefaultPrintThresholds::diagnostics) {
    ExEnv::out0() << indent << "Collected contributions to the intermediates from all tasks";
    if (to_all_tasks)
      ExEnv::out0() << " and distributed to every task" << endl;
    else
      ExEnv::out0() << " on task 0" << endl;
  }
}

const Ref<MOIndexSpace>&
R12IntEval::occ_act(SpinCase1 S) const
{
  return r12info()->refinfo()->occ_act(S);
}

const Ref<MOIndexSpace>&
R12IntEval::occ(SpinCase1 S) const
{
  return r12info()->refinfo()->occ(S);
}

const Ref<MOIndexSpace>&
R12IntEval::vir_act(SpinCase1 S) const
{
  if (r12info()->basis_vir() != r12info()->refinfo()->ref()->basis())
    return r12info()->vir_act(S);
  else
    return r12info()->refinfo()->uocc_act(S);
}

const Ref<MOIndexSpace>&
R12IntEval::vir(SpinCase1 S) const
{
  if (r12info()->basis_vir() != r12info()->refinfo()->ref()->basis())
    return r12info()->vir(S);
  else
    return r12info()->refinfo()->uocc(S);
}

const Ref<MOIndexSpace>&
R12IntEval::xspace(SpinCase1 S) const
{
    switch(r12info()->ansatz()->orbital_product()) {
    case LinearR12::OrbProd_ij:
	return(occ_act(S));
    case LinearR12::OrbProd_pq:
	return(r12info()->refinfo()->orbs(S));
    default:
	throw ProgrammingError("R12IntEval::xspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}

const Ref<MOIndexSpace>&
R12IntEval::hj_x_P(SpinCase1 S)
{
    switch(r12info()->ansatz()->orbital_product()) {
    case LinearR12::OrbProd_ij:
	return(hj_i_P(S));
    case LinearR12::OrbProd_pq:
	return(hj_p_P(S));
    default:
	throw ProgrammingError("R12IntEval::xspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}

const Ref<MOIndexSpace>&
R12IntEval::hj_x_p(SpinCase1 S)
{
    switch(r12info()->ansatz()->orbital_product()) {
    case LinearR12::OrbProd_ij:
	return(hj_i_p(S));
    case LinearR12::OrbProd_pq:
	return(hj_p_p(S));
    default:
	throw ProgrammingError("R12IntEval::xspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}

const Ref<MOIndexSpace>&
R12IntEval::hj_x_m(SpinCase1 S)
{
    switch(r12info()->ansatz()->orbital_product()) {
    case LinearR12::OrbProd_ij:
        return(hj_i_m(S));
    case LinearR12::OrbProd_pq:
        return(hj_p_m(S));
    default:
        throw ProgrammingError("R12IntEval::xspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}

const Ref<MOIndexSpace>&
R12IntEval::hj_x_a(SpinCase1 S)
{
    switch(r12info()->ansatz()->orbital_product()) {
    case LinearR12::OrbProd_ij:
        return(hj_i_a(S));
    case LinearR12::OrbProd_pq:
        return(hj_p_a(S));
    default:
        throw ProgrammingError("R12IntEval::xspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}

const Ref<MOIndexSpace>&
R12IntEval::hj_x_A(SpinCase1 S)
{
    switch(r12info()->ansatz()->orbital_product()) {
    case LinearR12::OrbProd_ij:
	return(hj_i_A(S));
    case LinearR12::OrbProd_pq:
	return(hj_p_A(S));
    default:
	throw ProgrammingError("R12IntEval::xspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}

const Ref<MOIndexSpace>&
R12IntEval::K_x_P(SpinCase1 S)
{
    switch(r12info()->ansatz()->orbital_product()) {
    case LinearR12::OrbProd_ij:
	return(K_i_P(S));
    case LinearR12::OrbProd_pq:
	return(K_p_P(S));
    default:
	throw ProgrammingError("R12IntEval::xspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}

const Ref<MOIndexSpace>&
R12IntEval::K_x_p(SpinCase1 S)
{
    switch(r12info()->ansatz()->orbital_product()) {
    case LinearR12::OrbProd_ij:
	return(K_i_p(S));
    case LinearR12::OrbProd_pq:
	return(K_p_p(S));
    default:
	throw ProgrammingError("R12IntEval::xspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}

const Ref<MOIndexSpace>&
R12IntEval::K_x_m(SpinCase1 S)
{
    switch(r12info()->ansatz()->orbital_product()) {
    case LinearR12::OrbProd_ij:
	return(K_i_m(S));
    case LinearR12::OrbProd_pq:
	return(K_p_m(S));
    default:
	throw ProgrammingError("R12IntEval::xspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}

const Ref<MOIndexSpace>&
R12IntEval::K_x_a(SpinCase1 S)
{
    switch(r12info()->ansatz()->orbital_product()) {
    case LinearR12::OrbProd_ij:
	return(K_i_a(S));
    case LinearR12::OrbProd_pq:
	return(K_p_a(S));
    default:
	throw ProgrammingError("R12IntEval::xspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}

const Ref<MOIndexSpace>&
R12IntEval::K_x_A(SpinCase1 S)
{
    switch(r12info()->ansatz()->orbital_product()) {
    case LinearR12::OrbProd_ij:
	return(K_i_A(S));
    case LinearR12::OrbProd_pq:
	return(K_p_A(S));
    default:
	throw ProgrammingError("R12IntEval::xspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}

const Ref<MOIndexSpace>&
R12IntEval::F_x_A(SpinCase1 S)
{
    switch(r12info()->ansatz()->orbital_product()) {
    case LinearR12::OrbProd_ij:
	return(F_i_A(S));
    case LinearR12::OrbProd_pq:
	return(F_p_A(S));
    default:
	throw ProgrammingError("R12IntEval::xspace() -- invalid orbital product of the R12 ansatz",__FILE__,__LINE__);
    }
}

namespace {
  /// Convert 2 spaces to SpinCase2
    SpinCase2
    spincase2(const Ref<MOIndexSpace>& space1,
              const Ref<MOIndexSpace>& space2)
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
      }
    }
};

std::string
R12IntEval::transform_label(const Ref<MOIndexSpace>& space1,
                            const Ref<MOIndexSpace>& space2,
                            const Ref<MOIndexSpace>& space3,
                            const Ref<MOIndexSpace>& space4) const
{
  std::ostringstream oss;
  // use physicists' notation
  oss << "<" << space1->id() << " " << space3->id() << "|" << space2->id() << " " << space4->id() << ">";
  // for case-insensitive file systems append spincase
  oss << "_" << id(spincase2(space1,space3));
  return oss.str();
}

std::string
R12IntEval::transform_label(const Ref<MOIndexSpace>& space1,
                            const Ref<MOIndexSpace>& space2,
                            const Ref<MOIndexSpace>& space3,
                            const Ref<MOIndexSpace>& space4,
                            unsigned int f12) const
{
  std::ostringstream oss;
  // use physicists' notation
  oss << "<" << space1->id() << " " << space3->id() << "| " << corrfactor()->label()
      << "[" << f12 << "] |" << space2->id() << " " << space4->id() << ">";
  // for case-insensitive file systems append spincase
  oss << "_" << id(spincase2(space1,space3));
  return oss.str();
}

std::string
R12IntEval::transform_label(const Ref<MOIndexSpace>& space1,
                            const Ref<MOIndexSpace>& space2,
                            const Ref<MOIndexSpace>& space3,
                            const Ref<MOIndexSpace>& space4,
                            unsigned int f12_left,
                            unsigned int f12_right) const
{
  std::ostringstream oss;
  // use physicists' notation
  oss << "<" << space1->id() << " " << space3->id() << "| " << corrfactor()->label()
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
