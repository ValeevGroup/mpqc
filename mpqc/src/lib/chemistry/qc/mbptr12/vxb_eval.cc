//
// vxb_eval.cc
//
// Copyright (C) 2003 Edward Valeev
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

#include <stdexcept>

#include <util/misc/formio.h>
#include <util/state/state_bin.h>
#include <util/ref/ref.h>
#include <math/scmat/local.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/vxb_eval.h>
#include <chemistry/qc/mbptr12/vxb_eval_sbs_a.h>
#include <chemistry/qc/mbptr12/vxb_eval_abs_a.h>
#include <chemistry/qc/mbptr12/vxb_eval_b.h>

using namespace std;
using namespace sc;

inline int max(int a,int b) { return (a > b) ? a : b;}

/*-----------
  R12IntEval
 -----------*/
static ClassDesc R12IntEval_cd(
  typeid(R12IntEval),"R12IntEval",1,"virtual public SavableState",
  0, 0, create<R12IntEval>);

R12IntEval::R12IntEval(MBPT2_R12* mbptr12)
{
  r12info_ = new R12IntEvalInfo(mbptr12);

  eval_sbs_a_ = new R12IntEval_sbs_A(r12info_);
  eval_abs_a_ = new R12IntEval_abs_A(r12info_);
  eval_b_ = new R12IntEval_B(r12info_);

  int nocc_act = r12info_->nocc_act();
  dim_aa_ = new SCDimension((nocc_act*(nocc_act-1))/2);
  dim_ab_ = new SCDimension(nocc_act*nocc_act);
  dim_s_ = new SCDimension((nocc_act*(nocc_act+1))/2);
  dim_t_ = dim_aa_;

  Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();
  Vaa_ = local_matrix_kit->matrix(dim_aa_,dim_aa_);
  Vab_ = local_matrix_kit->matrix(dim_ab_,dim_ab_);
  Xaa_ = local_matrix_kit->matrix(dim_aa_,dim_aa_);
  Xab_ = local_matrix_kit->matrix(dim_ab_,dim_ab_);
  Baa_ = local_matrix_kit->matrix(dim_aa_,dim_aa_);
  Bab_ = local_matrix_kit->matrix(dim_ab_,dim_ab_);
  emp2pair_aa_ = local_matrix_kit->vector(dim_aa_);
  emp2pair_ab_ = local_matrix_kit->vector(dim_ab_);

  gebc_ = mbptr12->gebc();
  stdapprox_ = mbptr12->stdapprox();
  spinadapted_ = mbptr12->spinadapted();

  // Default values
  evaluated_ = false;
  debug_ = 0;
}

R12IntEval::R12IntEval(StateIn& si) : SavableState(si)
{
  r12info_ << SavableState::restore_state(si);
  eval_sbs_a_ << SavableState::restore_state(si);
  eval_abs_a_ << SavableState::restore_state(si);
  eval_b_ << SavableState::restore_state(si);

  dim_aa_ << SavableState::restore_state(si);
  dim_ab_ << SavableState::restore_state(si);
  dim_s_ << SavableState::restore_state(si);
  dim_t_ << SavableState::restore_state(si);

  Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();
  Vaa_ = local_matrix_kit->matrix(dim_aa_,dim_aa_);
  Vab_ = local_matrix_kit->matrix(dim_ab_,dim_ab_);
  Xaa_ = local_matrix_kit->matrix(dim_aa_,dim_aa_);
  Xab_ = local_matrix_kit->matrix(dim_ab_,dim_ab_);
  Baa_ = local_matrix_kit->matrix(dim_aa_,dim_aa_);
  Bab_ = local_matrix_kit->matrix(dim_ab_,dim_ab_);
  emp2pair_aa_ = local_matrix_kit->vector(dim_aa_);
  emp2pair_ab_ = local_matrix_kit->vector(dim_ab_);

  Vaa_.restore(si);
  Vab_.restore(si);
  Xaa_.restore(si);
  Xab_.restore(si);
  Baa_.restore(si);
  Bab_.restore(si);
  emp2pair_aa_.restore(si);
  emp2pair_ab_.restore(si);

  int gebc; si.get(gebc); gebc_ = (bool) gebc;
  int stdapprox; si.get(stdapprox); stdapprox_ = (LinearR12::StandardApproximation) stdapprox;
  int spinadapted; si.get(spinadapted); spinadapted_ = (bool) spinadapted;
  int evaluated; si.get(evaluated); evaluated_ = (bool) evaluated;
  si.get(debug_);
}

R12IntEval::~R12IntEval()
{
  r12info_ = 0;
  eval_sbs_a_ = 0;
  eval_abs_a_ = 0;

  dim_aa_ = 0;
  dim_ab_ = 0;
  dim_s_ = 0;
  dim_t_ = 0;
}

void R12IntEval::save_data_state(StateOut& so)
{
  SavableState::save_state(r12info_.pointer(),so);
  SavableState::save_state(eval_sbs_a_.pointer(),so);
  SavableState::save_state(eval_abs_a_.pointer(),so);
  SavableState::save_state(eval_b_.pointer(),so);

  SavableState::save_state(dim_aa_.pointer(),so);
  SavableState::save_state(dim_ab_.pointer(),so);
  SavableState::save_state(dim_s_.pointer(),so);
  SavableState::save_state(dim_t_.pointer(),so);

  Vaa_.save(so);
  Vab_.save(so);
  Xaa_.save(so);
  Xab_.save(so);
  Baa_.save(so);
  Bab_.save(so);
  emp2pair_aa_.save(so);
  emp2pair_ab_.save(so);

  so.put((int)stdapprox_);
  so.put((int)spinadapted_);
  so.put((int)evaluated_);
  so.put(debug_);
}

void R12IntEval::obsolete()
{
  evaluated_ = false;
}

void R12IntEval::set_gebc(bool gebc) { gebc_ = gebc; };
void R12IntEval::set_stdapprox(LinearR12::StandardApproximation stdapprox) { stdapprox_ = stdapprox; };
void R12IntEval::set_spinadapted(bool spinadapted) { spinadapted_ = spinadapted; };
void R12IntEval::set_debug(int debug) { if (debug >= 0) { debug_ = debug; r12info_->set_debug_level(debug_); }};
void R12IntEval::set_dynamic(bool dynamic) { r12info_->set_dynamic(dynamic); };
void R12IntEval::set_print_percent(double pp) { r12info_->set_print_percent(pp); };
void R12IntEval::set_memory(size_t nbytes) { r12info_->set_memory(nbytes); };

Ref<R12IntEvalInfo> R12IntEval::r12info() const { return r12info_; };
RefSCDimension R12IntEval::dim_aa() const { return dim_aa_; };
RefSCDimension R12IntEval::dim_ab() const { return dim_ab_; };
RefSCDimension R12IntEval::dim_s() const { return dim_s_; };
RefSCDimension R12IntEval::dim_t() const { return dim_t_; };
RefDiagSCMatrix R12IntEval::evals() const { return r12info_->evals(); };

RefSCMatrix R12IntEval::V_aa() {
  compute();
  return Vaa_;
}

RefSCMatrix R12IntEval::X_aa() {
  compute();
  return Xaa_;
}

RefSCMatrix R12IntEval::B_aa() {
  compute();
  return Baa_;
}

RefSCMatrix R12IntEval::V_ab() {
  compute();
  return Vab_;
}

RefSCMatrix R12IntEval::X_ab() {
  compute();
  return Xab_;
}

RefSCMatrix R12IntEval::B_ab() {
  compute();
  return Bab_;
}

RefSCVector R12IntEval::emp2_aa() {
  compute();
  return emp2pair_aa_;
}

RefSCVector R12IntEval::emp2_ab() {
  compute();
  return emp2pair_ab_;
}

void R12IntEval::compute()
{
  if (evaluated_)
    return;
  
  if (gebc_)
    throw std::runtime_error("R12IntEval::compute: intermediates for MP2-R12 methods that assume neither GBC nor EBC have not been implemented yet");
  if (spinadapted_)
    throw std::runtime_error("R12IntEval::compute: spin-adapted R12 intermediates have not been implemented yet");

  int nocc_act = r12info()->nocc_act();
  int naa = nocc_act*(nocc_act-1)/2;
  int nab = nocc_act*nocc_act;
  int me = r12info()->msg()->me();
  MolecularEnergy* mole = r12info()->mole();

  if (me == 0 && mole->if_to_checkpoint()) {
    StateOutBin stateout(mole->checkpoint_file());
    SavableState::save_state(mole,stateout);
    ExEnv::out0() << indent << "Checkpointed the wave function" << endl;
  }

  // If doing standard MP2-R12/A or MP2-R12/A' methods contributions to the intermediates from different classes of
  // integrals can be computed separately ...
  if (gebc_ && stdapprox_ != LinearR12::StdApprox_B) {
    eval_sbs_a_->compute(Vaa_,Xaa_,Baa_,Vab_,Xab_,Bab_,emp2pair_aa_,emp2pair_ab_);
  
    if (r12info_->basis() != r12info_->basis_aux()) {
      if (me == 0 && debug_ > 1) {
        Vaa_.print("Alpha-alpha SBS V matrix");
        Baa_.print("Alpha-alpha SBS B matrix");
        Xaa_.print("Alpha-alpha SBS X matrix");
        Vab_.print("Alpha-beta SBS V matrix");
        Bab_.print("Alpha-beta SBS B matrix");
        Xab_.print("Alpha-beta SBS X matrix");
      }

      eval_abs_a_->compute(Vaa_,Xaa_,Baa_,Vab_,Xab_,Bab_);
    }
  }
  // ... otherwise (for MP2-R12/B or when BCs are not assumed) need to compute all classes
  // of integrals and then compute contributions in one shot. This is done by an R12IntEval_B object
  {
    eval_b_->compute(Vaa_,Xaa_,Baa_,Vab_,Xab_,Bab_,emp2pair_aa_,emp2pair_ab_);
  }

  if (me == 0 && debug_) {
    Vaa_.print("Alpha-alpha V matrix");
    Baa_.print("Alpha-alpha B matrix");
    Xaa_.print("Alpha-alpha X matrix");
    Vab_.print("Alpha-beta V matrix");
    Bab_.print("Alpha-beta B matrix");
    Xab_.print("Alpha-beta X matrix");
  }

  //
  // Compute basis set completeness
  //
  double traceV_aa = Vaa_->trace();
  double traceB_aa = Baa_->trace();
  double traceV_ab = Vab_->trace();
  double traceB_ab = Bab_->trace();
  
  ExEnv::out0() << endl;
  ExEnv::out0() << indent << "Basis Set completeness diagnostics:" << endl;
  ExEnv::out0() << indent << indent
		<< "-Tr(V)/Tr(B) for alpha-alpha pairs:" << indent <<
    scprintf("%10.6lf",(-1.0)*traceV_aa/traceB_aa) << endl;
  ExEnv::out0() << indent << indent
		<< "-Tr(V)/Tr(B) for alpha-beta pairs:" << indent <<
    scprintf("%10.6lf",(-1.0)*traceV_ab/traceB_ab) << endl;

  evaluated_ = true;

  if (me == 0 && mole->if_to_checkpoint()) {
    StateOutBin stateout(mole->checkpoint_file());
    SavableState::save_state(mole,stateout);
    ExEnv::out0() << indent << "Checkpointed the wave function" << endl;
  }
}
			 


/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
