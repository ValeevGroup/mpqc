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
#include <util/ref/ref.h>
#include <math/scmat/local.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/vxb_eval.h>
#include <chemistry/qc/mbptr12/vxb_eval_sbs_a.h>
#include <chemistry/qc/mbptr12/vxb_eval_abs_a.h>

using namespace std;
using namespace sc;

inline int max(int a,int b) { return (a > b) ? a : b;}

/*-----------
  R12IntEval
 -----------*/
static ClassDesc R12IntEval_cd(
  typeid(R12IntEval),"R12IntEval",1,"virtual public SavableState",
  0, 0, create<R12IntEval>);

R12IntEval::R12IntEval(MBPT2_R12& mbptr12)
{
  r12info_ = new R12IntEvalInfo(mbptr12);

  eval_sbs_a_ = new R12IntEval_sbs_A(r12info_);
  eval_abs_a_ = new R12IntEval_abs_A(r12info_);

  int nocc_act = r12info_->nocc_act();
  dim_aa_ = new SCDimension((nocc_act*(nocc_act-1))/2);
  dim_ab_ = new SCDimension(nocc_act*nocc_act);
  dim_s_ = new SCDimension((nocc_act*(nocc_act+1))/2);
  dim_t_ = dim_aa_;

  // Default values
  stdapprox_ = LinearR12::StdApprox_Ap;
  debug_ = 0;
}

R12IntEval::R12IntEval(StateIn& si) : SavableState(si)
{
  r12info_ << SavableState::restore_state(si);
  eval_sbs_a_ << SavableState::restore_state(si);
  eval_abs_a_ << SavableState::restore_state(si);

  int nocc_act = r12info_->nocc_act();
  dim_aa_ = new SCDimension((nocc_act*(nocc_act-1))/2);
  dim_ab_ = new SCDimension(nocc_act*nocc_act);
  dim_s_ = new SCDimension((nocc_act*(nocc_act+1))/2);
  dim_t_ = dim_aa_;

  int stdapprox;
  si.get(stdapprox);
  stdapprox_ = (LinearR12::StandardApproximation) stdapprox;
  spinadapted_ = false;
  debug_ = 0;
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

  so.put(stdapprox_);
}


void R12IntEval::compute(RefSCMatrix& Vaa,
			 RefSCMatrix& Xaa,
			 RefSCMatrix& Baa,
			 RefSCMatrix& Vab,
			 RefSCMatrix& Xab,
			 RefSCMatrix& Bab,
			 RefSCVector& emp2pair_aa,
			 RefSCVector& emp2pair_ab)
{
  if (spinadapted_)
    throw std::runtime_error("R12IntEval::compute: spin-adapted R12 intermediates have not been implemented yet");

  int nocc_act = r12info()->nocc_act();
  int naa = nocc_act*(nocc_act-1)/2;
  int nab = nocc_act*nocc_act;
  int me = r12info()->msg()->me();

  Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();
  if (me == 0) {
    Vaa = local_matrix_kit->matrix(dim_aa_,dim_aa_);
    Baa = local_matrix_kit->matrix(dim_aa_,dim_aa_);
    Xaa = local_matrix_kit->matrix(dim_aa_,dim_aa_);
    Vab = local_matrix_kit->matrix(dim_ab_,dim_ab_);
    Bab = local_matrix_kit->matrix(dim_ab_,dim_ab_);
    Xab = local_matrix_kit->matrix(dim_ab_,dim_ab_);
    emp2pair_aa = local_matrix_kit->vector(dim_aa_);
    emp2pair_ab = local_matrix_kit->vector(dim_ab_);
    Vaa->unit();
    Vab->unit();
    Baa->unit();
    Bab->unit();
    Xaa.assign(0.0);
    Xab.assign(0.0);
    emp2pair_aa.assign(0.0);
    emp2pair_ab.assign(0.0);

    if (debug_)
      ExEnv::out0() << indent << "Allocated V, B, and X intermediates" << endl;
    
    eval_sbs_a_->compute(Vaa,Xaa,Baa,Vab,Xab,Bab,emp2pair_aa,emp2pair_ab);
    eval_abs_a_->compute(Vaa,Xaa,Baa,Vab,Xab,Bab);

    if (debug_) {
      Vaa.print("Alpha-alpha V matrix");
      Baa.print("Alpha-alpha B matrix");
      Xaa.print("Alpha-alpha X matrix");
      Vab.print("Alpha-beta V matrix");
      Bab.print("Alpha-beta B matrix");
      Xab.print("Alpha-beta X matrix");
    }

    //
    // Compute basis set completeness
    //
    double traceV_aa = Vaa->trace();
    double traceB_aa = Baa->trace();
    double traceV_ab = Vab->trace();
    double traceB_ab = Bab->trace();

    ExEnv::out0() << endl;
    ExEnv::out0() << indent << "Basis Set completeness diagnostics:" << endl;
    ExEnv::out0() << indent
		 << "-Tr(V)/Tr(B) for alpha-alpha pairs:" << indent <<
      scprintf("%10.6lf",(-1.0)*traceV_aa/traceB_aa) << endl;
    ExEnv::out0() << indent
		 << "-Tr(V)/Tr(B) for alpha-beta pairs:" << indent <<
      scprintf("%10.6lf",(-1.0)*traceV_ab/traceB_ab) << endl;

  }
    
}
			 


/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
