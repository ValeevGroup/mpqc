//
// compute_energy_a.cc
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

#include <stdexcept>
#include <scconfig.h>
#include <util/misc/math.h>
#include <util/misc/formio.h>
#include <util/misc/timer.h>
#include <math/scmat/abstract.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/mp2r12_energy.h>
#include <chemistry/qc/mbptr12/transform_factory.h>

using namespace std;
using namespace sc;

void
MBPT2_R12::compute_energy_()
{
  tim_enter("mp2-f12 energy");

  //int DebugWait = 1;
  //while (DebugWait) {}
  
  Ref<R12IntEvalInfo> r12info;
  if (r12eval_.null()) {
    r12info = new R12IntEvalInfo(this);
    r12info->set_dynamic(dynamic_);
    r12info->set_print_percent(print_percent_);
    r12info->set_memory(mem_alloc);
    r12eval_ = new R12IntEval(r12info);
    r12eval_->set_debug(debug_);
  }
  else
    r12info = r12eval_->r12info();
  // This will actually compute the intermediates
  r12eval_->compute();
  
  double etotal = 0.0;
  
  //
  // Now we can compute and print pair energies
  //
  
  // can use projector 3 only for app C
  if (r12info->ansatz()->projector() != LinearR12::Projector_3) {

    // MP2-F12/A
    if (stdapprox_ == LinearR12::StdApprox_A ||
        stdapprox_ == LinearR12::StdApprox_Ap ||
        stdapprox_ == LinearR12::StdApprox_B) {
      tim_enter("mp2-f12/a pair energies");
      if (r12a_energy_.null())
        r12a_energy_ = new MP2R12Energy(r12eval_,LinearR12::StdApprox_A,debug_);
      r12a_energy_->print_pair_energies(spinadapted_);
      etotal = r12a_energy_->energy();
      tim_exit("mp2-f12/a pair energies");
    }
    
    // MP2-F12/A'
    // skip if using diagonal ansatz -- A' is equivalent to A then
    const bool skip_Ap = r12info->ansatz()->diag();
    if (!skip_Ap &&
	stdapprox_ == LinearR12::StdApprox_Ap ||
        stdapprox_ == LinearR12::StdApprox_B) {
      tim_enter("mp2-f12/a' pair energies");
      if (r12ap_energy_.null())
        r12ap_energy_ = new MP2R12Energy(r12eval_,LinearR12::StdApprox_Ap,debug_);
      r12ap_energy_->print_pair_energies(spinadapted_);
      etotal = r12ap_energy_->energy();
      tim_exit("mp2-f12/a' pair energies");
    }

    // MP2-F12/B
    if (stdapprox_ == LinearR12::StdApprox_B) {
      tim_enter("mp2-f12/b pair energies");
      if (r12b_energy_.null())
        r12b_energy_ = new MP2R12Energy(r12eval_,LinearR12::StdApprox_B,debug_);
      r12b_energy_->print_pair_energies(spinadapted_);
      etotal = r12b_energy_->energy();
      tim_exit("mp2-f12/b pair energies");
    }
    
    // MP2-F12/A''
    if (stdapprox_ == LinearR12::StdApprox_App) {
      tim_enter("mp2-f12/a'' pair energies");
      if (r12app_energy_.null())
        r12app_energy_ = new MP2R12Energy(r12eval_,LinearR12::StdApprox_App,debug_);
      r12app_energy_->print_pair_energies(spinadapted_);
      etotal = r12app_energy_->energy();
      tim_exit("mp2-f12/a'' pair energies");
    }

  } // end of != ansatz_3
  
  // MP2-F12/C
  if (stdapprox_ == LinearR12::StdApprox_C) {
    tim_enter("mp2-f12/c pair energies");
    if (r12c_energy_.null())
      r12c_energy_ = new MP2R12Energy(r12eval_,LinearR12::StdApprox_C,debug_);
    r12c_energy_->print_pair_energies(spinadapted_);
    etotal = r12c_energy_->energy();
    tim_exit("mp2-f12/c pair energies");
  }
  
  tim_exit("mp2-f12 energy");

  etotal += ref_energy();
  set_energy(etotal);
  set_actual_value_accuracy(reference_->actual_value_accuracy()
                            *ref_to_mp2_acc);
  
#if MP2R12ENERGY_CAN_COMPUTE_PAIRFUNCTION
  if (twopdm_grid_.nonnull()) {
    Ref<MP2R12Energy> wfn_to_plot;
    switch(stdapprox()) {
      case LinearR12::StdApprox_A:    wfn_to_plot = r12a_energy_;   break;
      case LinearR12::StdApprox_Ap:   wfn_to_plot = r12ap_energy_;  break;
      case LinearR12::StdApprox_App:  wfn_to_plot = r12app_energy_; break;
      case LinearR12::StdApprox_B:    wfn_to_plot = r12b_energy_;   break;
      case LinearR12::StdApprox_C:    wfn_to_plot = r12c_energy_;   break;
    }
    for(unsigned int sc2=0; sc2<NSpinCases2; sc2++) {
      wfn_to_plot->compute_pair_function(plot_pair_function_[0],plot_pair_function_[1],
                                         static_cast<SpinCase2>(sc2),
                                         twopdm_grid_);
    }
  }
#endif
    
  
  return;
}


////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
