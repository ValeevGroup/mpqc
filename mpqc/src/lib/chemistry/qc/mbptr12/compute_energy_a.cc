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
#include <util/misc/formio.h>
#include <util/misc/timer.h>
#include <math/scmat/abstract.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/vxb_eval.h>

using namespace std;
using namespace sc;


void
MBPT2_R12::compute_energy_a_spinorb_()
{
  tim_enter("mp2-r12a energy");

  double escf = ref_energy();

  int nocc = 0;
  for (int i=0; i<oso_dimension()->n(); i++) {
    if (reference_->occupation(i) == 2.0) nocc++;
  }
  int nocc_act = nocc - nfzc;
  int me = msg_->me();
  if (me == 0) {

    Ref<R12IntEval> r12eval = new R12IntEval(*this);
    r12eval->set_stdapprox(stdapprox_);
    // will spin-adapt the energies rather than the intermediates
    r12eval->set_spinadapted(false);
    r12eval->set_ints_file(r12ints_file_);
    r12eval->set_debug(debug_);
    r12eval->set_dynamic(dynamic_);
    r12eval->set_memory(mem_alloc);

    RefSCMatrix Vaa, Xaa, Baa, Vab, Xab, Bab;
    RefSCVector emp2_aa, emp2_ab;
    r12eval->compute(Vaa,Xaa,Baa,Vab,Xab,Bab,emp2_aa,emp2_ab);

    // Need eigenvalues
    RefDiagSCMatrix evalmat = r12eval->evals();
    double* evals = new double[nocc_act];
    for(int i=nfzc; i<nocc; i++)
      evals[i-nfzc] = evalmat(i);
    evalmat = 0;

    //
    // Evaluate pair energies
    //
    Ref<SCMatrixKit> localkit = Vaa.kit();
    RefSCDimension dim_aa = r12eval->dim_aa();
    RefSCDimension dim_ab = r12eval->dim_ab();
    RefSCVector epair_aa = localkit->vector(dim_aa);
    RefSCVector epair_ab = localkit->vector(dim_ab);
    RefSCVector er12_aa = localkit->vector(dim_aa);
    RefSCVector er12_ab = localkit->vector(dim_ab);
    double emp2tot_aa = 0.0;
    double emp2tot_ab = 0.0;
    double er12tot_aa = 0.0;
    double er12tot_ab = 0.0;
    
    // Alpha-alpha pairs
    RefSCMatrix Baa_ij = localkit->matrix(dim_aa,dim_aa);
    int ij=0;
    for(int i=0; i<nocc_act; i++)
      for(int j=0; j<i; j++, ij++) {

      RefSCVector Vaa_ij = Vaa.get_column(ij);

      // Form B(ij)kl,ow = Bkl,ow + 1/2(ek + el + eo + ew - 2ei - 2ej)Xkl,ow
      Baa_ij.assign(Baa);
      if (stdapprox_ == LinearR12::StdApprox_Ap) { 
      int kl=0;
      for(int k=0; k<nocc_act; k++)
	for(int l=0; l<k; l++, kl++) {
	  int ow=0;
	  for(int o=0; o<nocc_act; o++)
	    for(int w=0; w<o; w++, ow++) {
	      double fx = 0.5 * (evals[k] + evals[l] + evals[o] + evals[w] - 2.0*evals[i] - 2.0*evals[j]) *
		Xaa.get_element(kl,ow);
	      Baa_ij.accumulate_element(kl,ow,fx);
	    }
	}
      if (debug_ > 1) {
	Baa_ij.print("Full A' alpha-alpha B matrix");
      }
      }
      // For some reason invert_this doesn't work here
      Baa_ij->gen_invert_this();

      if (debug_ > 1) {
	Baa_ij.print("Inverse alpha-alpha B matrix");
      }

      double eaa_ij = -2.0*Vaa_ij.dot(Baa_ij * Vaa_ij);
      er12tot_aa += eaa_ij;
      er12_aa.set_element(ij,eaa_ij);
      emp2tot_aa += emp2_aa.get_element(ij);
    }
    Baa_ij=0;
    epair_aa->assign(emp2_aa);
    epair_aa->accumulate(er12_aa);
  
    // Alpha-beta pairs
    RefSCMatrix Bab_ij = localkit->matrix(dim_ab,dim_ab);
    ij=0;
    for(int i=0; i<nocc_act; i++)
      for(int j=0; j<nocc_act; j++, ij++) {

      RefSCVector Vab_ij = Vab.get_column(ij);

      // Form B(ij)kl,ow = Bkl,ow + 1/2(ek + el + eo + ew - 2ei - 2ej)Xkl,ow
      Bab_ij.assign(Bab);
      if (stdapprox_ == LinearR12::StdApprox_Ap) {
      int kl=0;
      for(int k=0; k<nocc_act; k++)
	for(int l=0; l<nocc_act; l++, kl++) {
	  int ow=0;
	  for(int o=0; o<nocc_act; o++)
	    for(int w=0; w<nocc_act; w++, ow++) {
	      double fx = 0.5 * (evals[k] + evals[l] + evals[o] + evals[w] - 2.0*evals[i] - 2.0*evals[j]) *
		Xab.get_element(kl,ow);
	      Bab_ij.accumulate_element(kl,ow,fx);
	    }
	}
      if (debug_ > 1) {
	Bab_ij.print("Full A' alpha-beta B matrix");
      }
      }
      // For some reason invert_this doesn't work here
      Bab_ij->gen_invert_this();

      if (debug_ > 1) {
	Bab_ij.print("Inverse alpha-beta B matrix");
      }

      double eab_ij = -1.0*Vab_ij.dot(Bab_ij * Vab_ij);
      er12tot_ab += eab_ij;
      er12_ab.set_element(ij,eab_ij);
      emp2tot_ab += emp2_ab.get_element(ij);
    }
    Bab_ij=0;
    epair_ab->assign(emp2_ab);
    epair_ab->accumulate(er12_ab);


    /*---------------------------------------
      Spin-adapt pair energies, if necessary
     ---------------------------------------*/
    if (!spinadapted_) {
      epair_0_ = epair_ab;
      epair_1_ = epair_aa;
    }
    else {
      epair_0_ = localkit->vector(r12eval->dim_s());
      epair_1_ = localkit->vector(r12eval->dim_t());

      // Triplet pairs are easy
      epair_1_->assign(epair_aa);
      epair_1_->scale(1.5);

      // Singlet pairs are trickier
      int ij_s = 0;
      for(int i=0; i<nocc_act; i++)
	for(int j=0; j<=i; j++, ij_s++) {
	  int ij_ab = i*nocc_act + j;
	  double eab = epair_ab->get_element(ij_ab);
	  epair_0_->set_element(ij_s,eab);
	}
    }

    if (stdapprox_ == LinearR12::StdApprox_A)
      ExEnv::out0() << endl << indent << "Alpha-alpha MBPT2-R12/A pair energies:" << endl;
    else
      ExEnv::out0() << endl << indent << "Alpha-alpha MBPT2-R12/A' pair energies:" << endl;
    ExEnv::out0() << indent << scprintf("    i       j        mp2(ij)        r12(ij)      mp2-r12(ij)") << endl;
    ExEnv::out0() << indent << scprintf("  -----   -----   ------------   ------------   ------------") << endl;
    for(int i=0,ij=0;i<nocc_act;i++)
      for(int j=0;j<i;j++,ij++) {
	ExEnv::out0() << indent << scprintf("  %3d     %3d     %12.9lf   %12.9lf   %12.9lf",i+1,j+1,
					    emp2_aa->get_element(ij),
					    er12_aa->get_element(ij),
					    epair_aa->get_element(ij)) << endl;
      }
    
    if (stdapprox_ == LinearR12::StdApprox_A)
      ExEnv::out0() << endl << indent << "Alpha-beta MBPT2-R12/A pair energies:" << endl;
    else
      ExEnv::out0() << endl << indent << "Alpha-beta MBPT2-R12/A' pair energies:" << endl;
    ExEnv::out0() << indent << scprintf("    i       j        mp2(ij)        r12(ij)      mp2-r12(ij)") << endl;
    ExEnv::out0() << indent << scprintf("  -----   -----   ------------   ------------   ------------") << endl;
    for(int i=0,ij=0;i<nocc_act;i++)
      for(int j=0;j<nocc_act;j++,ij++) {
	ExEnv::out0() << indent << scprintf("  %3d     %3d     %12.9lf   %12.9lf   %12.9lf",i+1,j+1,
					    emp2_ab->get_element(ij),
					    er12_ab->get_element(ij),
					    epair_ab->get_element(ij)) << endl;
      }

    mp2_corr_energy_ = emp2tot_aa + emp2tot_ab;
    r12_corr_energy_ = er12tot_aa + er12tot_ab;

    delete[] evals;
  }
  tim_exit("mp2-r12a energy");
  
  ///////////////////////////////////////////////////////////////
  // The computation of the MP2 energy is now complete on each
  // node;
  ///////////////////////////////////////////////////////////////

  double etotal = escf + mp2_corr_energy_ + r12_corr_energy_;
  ExEnv::out0()<<endl<<indent
	       <<scprintf("RHF energy [au]:                           %17.12lf\n", escf);
  ExEnv::out0()<<indent
	       <<scprintf("MP2 correlation energy [au]:               %17.12lf\n", mp2_corr_energy_);
  if (stdapprox_ == LinearR12::StdApprox_A)
    ExEnv::out0()<<indent
		 <<scprintf("(MBPT2)-R12/A  correlation energy [au]:    %17.12lf\n", r12_corr_energy_);
  else
    ExEnv::out0()<<indent
		 <<scprintf("(MBPT2)-R12/A' correlation energy [au]:    %17.12lf\n", r12_corr_energy_);
  if (stdapprox_ == LinearR12::StdApprox_A)
    ExEnv::out0()<<indent
		 <<scprintf("MBPT2-R12/A  correlation energy [au]:      %17.12lf\n", mp2_corr_energy_ + r12_corr_energy_);
  else
    ExEnv::out0()<<indent
		 <<scprintf("MBPT2-R12/A' correlation energy [au]:      %17.12lf\n", mp2_corr_energy_ + r12_corr_energy_);
  if (stdapprox_ == LinearR12::StdApprox_A)
    ExEnv::out0()<<indent
		 <<scprintf("MBPT2-R12/A  energy [au]:                  %17.12lf\n", etotal);
  else
    ExEnv::out0()<<indent
		 <<scprintf("MBPT2-R12/A' energy [au]:                  %17.12lf\n", etotal);
  ExEnv::out0().flush();

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
