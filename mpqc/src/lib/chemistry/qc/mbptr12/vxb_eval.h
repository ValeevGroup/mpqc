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

class MBPT2_R12;
class R12IntEvalInfo;
class R12IntEval_sbs_A;
class R12IntEval_abs_A;

  /** Class R12IntEval is the top-level R12 intermediate evaluator */

class R12IntEval : virtual public SavableState {

  /// This describes the state of the evaluator - whether it's been evaluated or not
  bool evaluated_;
  
  Ref<R12IntEvalInfo> r12info_;

  Ref<R12IntEval_sbs_A> eval_sbs_a_;
  Ref<R12IntEval_abs_A> eval_abs_a_;

  RefSCMatrix Vaa_, Vab_, Xaa_, Xab_, Baa_, Bab_;
  RefSCVector emp2pair_aa_, emp2pair_ab_;
  RefSCDimension dim_aa_, dim_ab_, dim_s_, dim_t_;

  LinearR12::StandardApproximation stdapprox_;
  bool spinadapted_;
  int debug_;

public:

  R12IntEval(StateIn&);
  R12IntEval(MBPT2_R12*);
  ~R12IntEval();

  void save_data_state(StateOut&);
  /// Makes the evaluator obsolete, the next call to compute() will cause the intermediates to be recomputed
  void obsolete();
  
  void set_stdapprox(LinearR12::StandardApproximation stdapprox);
  void set_spinadapted(bool spinadapted);
  void set_debug(int debug);
  void set_dynamic(bool dynamic);
  void set_print_percent(double print_percent);
  void set_memory(size_t nbytes);

  Ref<R12IntEvalInfo> r12info() const;
  RefSCDimension dim_aa() const;
  RefSCDimension dim_ab() const;
  RefSCDimension dim_s() const;
  RefSCDimension dim_t() const;

  RefSCMatrix V_aa();
  RefSCMatrix X_aa();
  RefSCMatrix B_aa();
  RefSCMatrix V_ab();
  RefSCMatrix X_ab();
  RefSCMatrix B_ab();
  RefSCVector emp2_aa();
  RefSCVector emp2_ab();
  
  void compute();
  RefDiagSCMatrix evals() const;
	       
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:


