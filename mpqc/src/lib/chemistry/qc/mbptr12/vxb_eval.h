//
// vxb_eval.h
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
#pragma interface
#endif

#ifndef _chemistry_qc_mbptr12_vxbeval_h
#define _chemistry_qc_mbptr12_vxbeval_h

#include <util/ref/ref.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/linearr12.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>
#include <chemistry/qc/mbptr12/vxb_eval_sbs_a.h>
#include <chemistry/qc/mbptr12/vxb_eval_abs_a.h>

namespace sc {


  /** Class R12IntEval is the top-level R12 intermediate evaluator */

class R12IntEval : virtual public SavableState {

  Ref<R12IntEvalInfo> r12info_;

  Ref<R12IntEval_sbs_A> eval_sbs_a_;
  Ref<R12IntEval_abs_A> eval_abs_a_;

  RefSCDimension dim_aa_, dim_ab_, dim_s_, dim_t_;

  LinearR12::StandardApproximation stdapprox_;
  bool spinadapted_;
  int debug_;

public:

  R12IntEval(StateIn&);
  R12IntEval(MBPT2_R12&);
  ~R12IntEval();

  void save_data_state(StateOut&);

  void set_stdapprox(LinearR12::StandardApproximation stdapprox) { stdapprox_ = stdapprox; };
  void set_spinadapted(bool spinadapted) { spinadapted_ = spinadapted; };
  void set_debug(int debug) { if (debug >= 0) { debug_ = debug; r12info_->set_debug_level(debug_); }};
  void set_dynamic(bool dynamic) { r12info_->set_dynamic(dynamic); };
  void set_checkpoint(bool chkpt) { r12info_->set_checkpoint(chkpt); };
  void set_checkpoint_file(const char* filename) { r12info_->set_checkpoint_file(filename); };
  void set_memory(size_t nbytes) { r12info_->set_memory(nbytes); };
  void set_ints_file(const char* filename) { r12info_->set_ints_file(filename); };

  Ref<R12IntEvalInfo> r12info() const { return r12info_; };
  RefSCDimension dim_aa() const { return dim_aa_; };
  RefSCDimension dim_ab() const { return dim_ab_; };
  RefSCDimension dim_s() const { return dim_s_; };
  RefSCDimension dim_t() const { return dim_t_; };

  void compute(RefSCMatrix& Vaa,
	       RefSCMatrix& Xaa,
	       RefSCMatrix& Baa,
	       RefSCMatrix& Vab,
	       RefSCMatrix& Xab,
	       RefSCMatrix& Bab,
	       RefSCVector& emp2pair_aa,
	       RefSCVector& emp2pair_ab);
  RefDiagSCMatrix evals() const { return r12info_->evals(); };
	       
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:


