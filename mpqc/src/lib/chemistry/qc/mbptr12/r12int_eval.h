//
// r12int_eval.h
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
#pragma interface
#endif

#ifndef _chemistry_qc_mbptr12_r12inteval_h
#define _chemistry_qc_mbptr12_r12inteval_h

#include <util/ref/ref.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>
#include <chemistry/qc/mbptr12/linearr12.h>

namespace sc {

  /** R12IntEval is the top-level class which computes intermediates occuring in linear R12 theories.
      This class is used by all Wavefunction classes that implement linear R12 methods.
  */

class R12IntEval : virtual public SavableState {

  bool evaluated_;

  // Calculation information (number of basis functions, R12 approximation, etc.)
  Ref<R12IntEvalInfo> r12info_;

  RefSCMatrix Vaa_, Vab_, Xaa_, Xab_, Baa_, Bab_, Aaa_, Aab_;
  RefSCVector emp2pair_aa_, emp2pair_ab_;
  RefSCDimension dim_aa_, dim_ab_, dim_s_, dim_t_;

  bool gebc_;
  LinearR12::ABSMethod abs_method_;
  LinearR12::StandardApproximation stdapprox_;
  bool spinadapted_;
  int debug_;

  // (act.occ. OBS| act.occ. OBS) transform
  Ref<TwoBodyMOIntsTransform> ipjq_tform_;
  // (act.occ. occ.| act.occ. RI-BS) transform
  Ref<TwoBodyMOIntsTransform> ikjy_tform_;

  /// Initialize transforms
  void init_tforms_();
  /// Set intermediates to zero + add the "diagonal" contributions
  void init_intermeds_();
  /// Compute r^2 contribution to X
  void r2_contrib_to_X_();  
  /// Checkpoint the top-level molecular energy
  void checkpoint_() const;

  /** Compute the number tasks which have access to the integrals.
      map_to_twi is a vector<int> each element of which will hold the
      number of tasks with access to the integrals of lower rank than that task
      (or -1 if the task doesn't have access to the integrals) */
  const int tasks_with_ints_(const Ref<R12IntsAcc> ints_acc, vector<int>& map_to_twi);
  
  /// Compute OBS contribution to V, X, and B (these contributions are independent of the method)
  void obs_contrib_to_VXB_gebc_();

  /// Compute 1-ABS contribution to V, X, and B (these contributions are independent of the method)
  void abs1_contrib_to_VXB_gebc_();
  
  /** Sum contributions to the intermediates from all nodes and broadcast so
      every node has the correct matrices */
  void globally_sum_intermeds_();
  
public:
  R12IntEval(StateIn&);
  R12IntEval(const Ref<R12IntEvalInfo>&);
  ~R12IntEval();

  void save_data_state(StateOut&);
  virtual void obsolete();

  void set_gebc(bool gebc);
  void set_absmethod(LinearR12::ABSMethod abs_method);
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

  /// This function causes the intermediate matrices to be computed.
  virtual void compute();

  /// Returns alpha-alpha block of the V intermediate matrix.
  RefSCMatrix V_aa();
  /// Returns alpha-alpha block of the X intermediate matrix.
  RefSCMatrix X_aa();
  /// Returns alpha-alpha block of the B intermediate matrix.
  RefSCMatrix B_aa();
  /// Returns alpha-alpha block of the A intermediate matrix.
  RefSCMatrix A_aa();
  /// Returns alpha-beta block of the V intermediate matrix.
  RefSCMatrix V_ab();
  /// Returns alpha-beta block of the X intermediate matrix.
  RefSCMatrix X_ab();
  /// Returns alpha-beta block of the B intermediate matrix.
  RefSCMatrix B_ab();
  /// Returns alpha-beta block of the A intermediate matrix.
  RefSCMatrix A_ab();
  /// Returns alpha-alpha MP2 pair energies.
  RefSCVector emp2_aa();
  /// Returns alpha-beta MP2 pair energies.
  RefSCVector emp2_ab();

  RefDiagSCMatrix evals() const;  
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:


