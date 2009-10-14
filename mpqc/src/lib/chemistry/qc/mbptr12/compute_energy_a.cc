//
// compute_energy_a.cc
//
// Copyright (C) 2003 Edward Valeev
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

#include <stdexcept>
#include <scconfig.h>
#include <util/misc/math.h>
#include <util/misc/formio.h>
#include <util/misc/regtime.h>
#include <math/scmat/abstract.h>
#include <math/scmat/local.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/mp2r12_energy.h>
#include <chemistry/qc/mbptr12/transform_factory.h>

using namespace std;
using namespace sc;

namespace {
    // returns the sum of all elements of v
    double sum(const RefSCVector& v);
    // returns the total R12 energy correction
    double er12(const Ref<MP2R12Energy>& mp2r12_energy);
}

void
MBPT2_R12::compute_energy_()
{
  Timer tim("mp2-f12 energy");

  //int DebugWait = 1;
  //while (DebugWait) {}

  const Ref<R12IntEvalInfo>& r12info = r12evalinfo_;
  r12info->initialize();

  // test app. C general-density B builder by feeding HF densities to it
#define TEST_GENERAL_FOCK_BUILD 0
#if TEST_GENERAL_FOCK_BUILD
  {
    Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
    RefSymmSCMatrix P[2];
    for(int s=0; s<NSpinCases1; ++s) {
      const SpinCase1 spin = static_cast<SpinCase1>(s);
      Ref<OrbitalSpace> occ = occ(spin);
      Ref<OrbitalSpace> orbs = orbs(spin);
      std::vector<unsigned int> occ2orbs = (*orbs << *occ);
      P[s] = localkit->symmmatrix(new SCDimension(orbs->rank()));
      P[s].assign(0.0);
      for(std::vector<unsigned int>::iterator i=occ2orbs.begin();
          i!=occ2orbs.end();
          ++i)
        P[s](*i, *i) = 1.0;
    }
    r12info->set_opdm(P[0], P[1]);
  }
#endif

  if (r12eval_.null()) {
    // since r12intevalinfo uses this class' KeyVal to initialize, dynamic is set automatically
    r12evalinfo_->set_print_percent(print_percent_);
    r12evalinfo_->set_memory(mem_alloc);
    r12eval_ = new R12IntEval(r12evalinfo_);
    r12eval_->set_debug(debug_);
  }
  // This will actually compute the intermediates
  r12eval_->compute();

  Ref<R12EnergyIntermediates> r12intermediates;

  double etotal = 0.0;
  double ef12 = 0.0;

  //
  // Now we can compute and print pair energies
  //

  // can use projector 3 only for app C
  if (r12info->ansatz()->projector() != LinearR12::Projector_3) {

    // MP2-F12/A'
    if (r12info->stdapprox() == LinearR12::StdApprox_Ap ||
        r12info->stdapprox() == LinearR12::StdApprox_B) {
      Timer tim2("mp2-f12/a' pair energies");
      if (r12ap_energy_.null()){
        r12intermediates=new R12EnergyIntermediates(r12eval_,LinearR12::StdApprox_Ap);
        r12ap_energy_ = construct_MP2R12Energy(r12intermediates,debug_,new_energy_);
      }
      r12ap_energy_->print_pair_energies(r12info->spinadapted());
      etotal = r12ap_energy_->energy();
      ef12 = er12(r12ap_energy_);
    }

    // MP2-F12/B
    if (r12info->stdapprox() == LinearR12::StdApprox_B) {
      Timer tim2("mp2-f12/b pair energies");
      if (r12b_energy_.null()){
        r12intermediates=new R12EnergyIntermediates(r12eval_,LinearR12::StdApprox_B);
        r12b_energy_ = construct_MP2R12Energy(r12intermediates,debug_,new_energy_);
      }
      r12b_energy_->print_pair_energies(r12info->spinadapted());
      etotal = r12b_energy_->energy();
      ef12 = er12(r12b_energy_);
    }

    // MP2-F12/A''
    if (r12info->stdapprox() == LinearR12::StdApprox_App) {
      Timer tim2("mp2-f12/a'' pair energies");
      if (r12app_energy_.null()){
        r12intermediates=new R12EnergyIntermediates(r12eval_,LinearR12::StdApprox_App);
        r12app_energy_ = construct_MP2R12Energy(r12intermediates,debug_,new_energy_);
      }
      r12app_energy_->print_pair_energies(r12info->spinadapted());
      etotal = r12app_energy_->energy();
      ef12 = er12(r12app_energy_);
    }

  } // end of != ansatz_3

  // MP2-F12/C
  if (r12info->stdapprox() == LinearR12::StdApprox_C ||
      r12info->stdapprox() == LinearR12::StdApprox_Cp) {
    Timer tim2("mp2-f12/c pair energies");
    if (r12c_energy_.null()){
      r12intermediates=new R12EnergyIntermediates(r12eval_,r12info->stdapprox());
      r12c_energy_ = construct_MP2R12Energy(r12intermediates,debug_,new_energy_);
    }
    r12c_energy_->print_pair_energies(r12info->spinadapted());
    etotal = r12c_energy_->energy();
    ef12 = er12(r12c_energy_);
  }

  if (cabs_singles_ && (r12evalinfo_->obs_eq_ribs() == false)) {
    cabs_singles_energy_ = r12eval_->emp2_cabs_singles();
  }

  tim.exit();

  mp2_corr_energy_ = etotal - ef12;
  etotal += ref_energy();
  etotal += cabs_singles_energy_;
  set_energy(etotal);
  set_actual_value_accuracy(reference_->actual_value_accuracy()
                            *ref_to_mp2r12_acc());

#if MP2R12ENERGY_CAN_COMPUTE_PAIRFUNCTION
  if (twopdm_grid_.nonnull()) {
    Ref<MP2R12Energy> wfn_to_plot;
    switch(r12info->stdapprox()) {
      case LinearR12::StdApprox_Ap:   wfn_to_plot = r12ap_energy_;  break;
      case LinearR12::StdApprox_App:  wfn_to_plot = r12app_energy_; break;
      case LinearR12::StdApprox_B:    wfn_to_plot = r12b_energy_;   break;
      case LinearR12::StdApprox_C:    wfn_to_plot = r12c_energy_;   break;
      case LinearR12::StdApprox_Cp:   wfn_to_plot = r12c_energy_;   break;
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

namespace {
    double sum(const RefSCVector& v) {
	RefSCVector unit = v.clone();  unit.assign(1.0);
	return v.dot(unit);
    }
    double er12(const Ref<MP2R12Energy>& mp2r12_energy) {
	double result = 0.0;
	for(int spin=0; spin<NSpinCases2; spin++) {
	    const SpinCase2 spincase2 = static_cast<SpinCase2>(spin);
	    result += sum(mp2r12_energy->ef12(spincase2));
	}
	return result;
    }
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
