//
// vxb_eval_abs_a.h
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

#ifndef _chemistry_qc_mbptr12_vxbevalabsa_h
#define _chemistry_qc_mbptr12_vxbevalabsa_h

#include <util/ref/ref.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>

namespace sc {


  /** Class R12IntEval_abs_A evaluates R12 intermediates for the auxiliary-basis part of
      the two-basis MP2-R12 theory in standard approximations A and A' */

class R12IntEval_abs_A : virtual public SavableState {

  Ref<R12IntEvalInfo> r12info_;

  int current_orbital_;
  int restart_orbital_;

  /* utility functions */
  int compute_transform_batchsize_(size_t mem_alloc, size_t mem_static, int nocc_act, const int num_te_types);
  distsize_t compute_transform_dynamic_memory_(int ni, int nocc_act, const int num_te_types);

public:
  R12IntEval_abs_A(StateIn&);
  R12IntEval_abs_A(Ref<R12IntEvalInfo>&);
  ~R12IntEval_abs_A();

  void save_data_state(StateOut&);

  void compute(RefSCMatrix& Vaa,
	       RefSCMatrix& Xaa,
	       RefSCMatrix& Baa,
	       RefSCMatrix& Vab,
	       RefSCMatrix& Xab,
	       RefSCMatrix& Bab);

  Ref<R12IntEvalInfo> r12info() const { return r12info_; };
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:


