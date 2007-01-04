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
#include <util/misc/regtime.h>
#include <math/scmat/abstract.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/mp2r12_energy.h>

using namespace std;
using namespace sc;


void
MBPT2_R12::compute_energy_a_()
{
  Timer tim("mp2-r12/a energy");

  if (r12eval_.null()) {
    Ref<R12IntEvalInfo> r12info = new R12IntEvalInfo(this);
    r12info->set_dynamic(dynamic_);
    r12info->set_print_percent(print_percent_);
    r12info->set_memory(mem_alloc);
    r12eval_ = new R12IntEval(r12info,gbc_,ebc_,abs_method_,stdapprox_);
    r12eval_->include_mp1(include_mp1_);
    r12eval_->set_debug(debug_);
  }
  // This will actually compute the intermediates
  r12eval_->compute();

  double etotal = 0.0;
  
  // Now we can compute and print pair energies
  tim.enter("mp2-r12/a pair energies");
  if (r12a_energy_.null())
    r12a_energy_ = new MP2R12Energy(r12eval_,LinearR12::StdApprox_A,debug_);
  r12a_energy_->print_pair_energies(spinadapted_);
  etotal = r12a_energy_->energy();
  tim.exit("mp2-r12/a pair energies");
  if (stdapprox_ == LinearR12::StdApprox_Ap) {
    tim.enter("mp2-r12/a' pair energies");
    if (r12ap_energy_.null())
      r12ap_energy_ = new MP2R12Energy(r12eval_,LinearR12::StdApprox_Ap,debug_);
    r12ap_energy_->print_pair_energies(spinadapted_);

    /*const double radius = 1.0;
    SCVector3 r1(0.0,0.0,radius);
    r12ap_energy_->compute_pair_function_ab(0,r1,r1);
    ExEnv::out0() << endl<<endl;
    const int nintervals = 100;
    const double Phi_start = -(M_PI);
    const double Phi_end = M_PI;
    const double dPhi = (Phi_end - Phi_start) / nintervals;
    const int npts = nintervals + 1;
    for(int i=-nintervals/2; i<=nintervals/2; i++) {
      const double Phi = i*dPhi;
      const double z = radius * cos(Phi);
      const double x = radius * sin(Phi);
      SCVector3 r2(x,0.0,z);
      ExEnv::out0() << indent << Phi;
      r12ap_energy_->compute_pair_function_ab(0,r1,r2);
      }*/
    if (twopdm_grid_aa_.nonnull())
      r12ap_energy_->compute_pair_function_aa(0,twopdm_grid_aa_);
    if (twopdm_grid_ab_.nonnull())
      r12ap_energy_->compute_pair_function_ab(0,twopdm_grid_ab_);

    etotal = r12ap_energy_->energy();
    tim.exit("mp2-r12/a' pair energies");
  }

  tim.exit("mp2-r12/a energy");

  etotal += ref_energy();
  set_energy(etotal);
  set_actual_value_accuracy(reference_->actual_value_accuracy()
                            *ref_to_mp2_acc);
  
  return;
}


////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
