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
#include <chemistry/qc/mbpt/bzerofast.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/vxb_eval.h>

using namespace std;
using namespace sc;


void
MBPT2_R12::compute_energy_a_()
{
  tim_enter("mp2-r12a energy");

  char * SA_str = (stdapprox_ == LinearR12::StdApprox_A) ? strdup("A") : strdup("A'");

  double escf = ref_energy();
  double emp2tot_aa = 0.0;
  double emp2tot_ab = 0.0;
  double er12tot_aa = 0.0;
  double er12tot_ab = 0.0;
  double emp2tot_0 = 0.0;
  double emp2tot_1 = 0.0;
  double er12tot_0 = 0.0;
  double er12tot_1 = 0.0;

  int nocc = 0;
  for (int i=0; i<oso_dimension()->n(); i++) {
    if (reference_->occupation(i) == 2.0) nocc++;
  }
  int nocc_act = nocc - nfzc;
  int me = msg_->me();
  int ntasks = msg_->n();

  r12eval_ = new R12IntEval(this);
  r12eval_->set_stdapprox(stdapprox_);
  // will spin-adapt the energies rather than the intermediates
  r12eval_->set_spinadapted(false);
  r12eval_->set_ints_file(r12ints_file_);
  r12eval_->set_debug(debug_);
  r12eval_->set_dynamic(dynamic_);
  r12eval_->set_memory(mem_alloc);
  
  RefSCMatrix Vaa, Xaa, Baa, Vab, Xab, Bab;
  RefSCVector emp2_aa, emp2_ab;
  r12eval_->compute(Vaa,Xaa,Baa,Vab,Xab,Bab,emp2_aa,emp2_ab);

  //
  // Evaluate pair energies:
  // distribute workload among nodes by pair index
  //

  // Need eigenvalues
  RefDiagSCMatrix evalmat = r12eval_->evals();
  double* evals = new double[nocc_act];
  for(int i=nfzc; i<nocc; i++)
    evals[i-nfzc] = evalmat(i);
  evalmat = 0;

  Ref<SCMatrixKit> localkit = Vaa.kit();
  RefSCDimension dim_aa = r12eval_->dim_aa();
  RefSCDimension dim_ab = r12eval_->dim_ab();
  int naa = dim_aa.n();
  int nab = dim_ab.n();
  RefSCVector epair_aa = localkit->vector(dim_aa);
  RefSCVector epair_ab = localkit->vector(dim_ab);
  RefSCVector er12_aa = localkit->vector(dim_aa);
  RefSCVector er12_ab = localkit->vector(dim_ab);
  double* er12_aa_vec = new double[naa];
  double* er12_ab_vec = new double[nab];

  // Alpha-alpha pairs
  RefSCMatrix Baa_ij = localkit->matrix(dim_aa,dim_aa);

  // In MP2-R12/A the B matrix is the same for all pairs
  if (stdapprox_ == LinearR12::StdApprox_A) {
    Baa_ij.assign(Baa);
    if (debug_ > 1)
      Baa_ij.print("Inverse alpha-alpha B matrix");
    Baa_ij->gen_invert_this();
  }

  bzerofast(er12_aa_vec,naa);
  int ij=0;
  for(int i=0; i<nocc_act; i++)
    for(int j=0; j<i; j++, ij++) {

      if (ij%ntasks != me)
        continue;

      RefSCVector Vaa_ij = Vaa.get_column(ij);

      // In MP2-R12/A' matrices B are pair-specific:
      // Form B(ij)kl,ow = Bkl,ow + 1/2(ek + el + eo + ew - 2ei - 2ej)Xkl,ow
      if (stdapprox_ == LinearR12::StdApprox_Ap) {
        Baa_ij.assign(Baa);
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
        if (debug_ > 1)
          Baa_ij.print("Full A' alpha-alpha B matrix");

        Baa_ij->gen_invert_this();

        if (debug_ > 1)
          Baa_ij.print("Inverse A' alpha-alpha B matrix");
      }

      double eaa_ij = -2.0*Vaa_ij.dot(Baa_ij * Vaa_ij);
      er12_aa_vec[ij] = eaa_ij;
    }
  Baa_ij=0;
  msg_->sum(er12_aa_vec,naa,0,-1);
  er12_aa->assign(er12_aa_vec);
  epair_aa->assign(emp2_aa);
  epair_aa->accumulate(er12_aa);
  for(int ij=0; ij<naa; ij++) {
    er12tot_aa += er12_aa->get_element(ij);
    emp2tot_aa += emp2_aa->get_element(ij);
  }
  delete[] er12_aa_vec;

  // Alpha-beta pairs
  RefSCMatrix Bab_ij = localkit->matrix(dim_ab,dim_ab);

  // In MP2-R12/A the B matrix is the same for all pairs
  if (stdapprox_ == LinearR12::StdApprox_A) {
    Bab_ij.assign(Bab);
    if (debug_ > 1)
      Bab_ij.print("Inverse alpha-beta B matrix");
    Bab_ij->gen_invert_this();
  }

  bzerofast(er12_ab_vec,nab);
  ij=0;
  for(int i=0; i<nocc_act; i++)
    for(int j=0; j<nocc_act; j++, ij++) {

      if (ij%ntasks != me)
        continue;

      RefSCVector Vab_ij = Vab.get_column(ij);

      // In MP2-R12/A' matrices B are pair-specific:
      // Form B(ij)kl,ow = Bkl,ow + 1/2(ek + el + eo + ew - 2ei - 2ej)Xkl,ow
      if (stdapprox_ == LinearR12::StdApprox_Ap) {
        Bab_ij.assign(Bab);
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
        if (debug_ > 1)
	    Bab_ij.print("Full A' alpha-beta B matrix");

        Bab_ij->gen_invert_this();

        if (debug_ > 1)
          Bab_ij.print("Inverse A' alpha-beta B matrix");
      }

      double eab_ij = -1.0*Vab_ij.dot(Bab_ij * Vab_ij);
      er12_ab_vec[ij] = eab_ij;
    }
  Bab_ij=0;
  msg_->sum(er12_ab_vec,nab,0,-1);
  epair_ab->assign(emp2_ab);
  epair_ab->accumulate(er12_ab);
  er12_ab->assign(er12_ab_vec);
  for(int ij=0; ij<nab; ij++) {
    er12tot_ab += er12_ab->get_element(ij);
    emp2tot_ab += emp2_ab->get_element(ij);
  }
  delete[] er12_ab_vec;

  /*---------------------------------------
    Spin-adapt pair energies, if necessary
   ---------------------------------------*/
  if (!spinadapted_) {
    epair_0_ = epair_ab;
    epair_1_ = epair_aa;

    ExEnv::out0() << endl << indent << "Alpha-alpha MBPT2-R12/" << SA_str << " pair energies:" << endl;
    ExEnv::out0() << indent << scprintf("    i       j        mp2(ij)        r12(ij)      mp2-r12(ij)") << endl;
    ExEnv::out0() << indent << scprintf("  -----   -----   ------------   ------------   ------------") << endl;
    for(int i=0,ij=0;i<nocc_act;i++)
      for(int j=0;j<i;j++,ij++) {
        ExEnv::out0() << indent << scprintf("  %3d     %3d     %12.9lf   %12.9lf   %12.9lf",i+1,j+1,
                                            emp2_aa->get_element(ij),
                                            er12_aa->get_element(ij),
                                            epair_aa->get_element(ij)) << endl;
      }
      
    ExEnv::out0() << endl << indent << "Alpha-beta MBPT2-R12/" << SA_str << " pair energies:" << endl;
    ExEnv::out0() << indent << scprintf("    i       j        mp2(ij)        r12(ij)      mp2-r12(ij)") << endl;
    ExEnv::out0() << indent << scprintf("  -----   -----   ------------   ------------   ------------") << endl;
    for(int i=0,ij=0;i<nocc_act;i++)
      for(int j=0;j<nocc_act;j++,ij++) {
        ExEnv::out0() << indent << scprintf("  %3d     %3d     %12.9lf   %12.9lf   %12.9lf",i+1,j+1,
                                            emp2_ab->get_element(ij),
                                            er12_ab->get_element(ij),
                                            epair_ab->get_element(ij)) << endl;
      }
  }
  else {
    epair_0_ = localkit->vector(r12eval_->dim_s());
    epair_1_ = localkit->vector(r12eval_->dim_t());
    RefSCVector emp2_0 = localkit->vector(r12eval_->dim_s());
    RefSCVector emp2_1 = localkit->vector(r12eval_->dim_t());
    RefSCVector er12_0 = localkit->vector(r12eval_->dim_s());
    RefSCVector er12_1 = localkit->vector(r12eval_->dim_t());

    // Triplet pairs are easy
    epair_1_->assign(epair_aa);
    epair_1_->scale(1.5);
    emp2_1->assign(emp2_aa);
    emp2_1->scale(1.5);
    er12_1->assign(er12_aa);
    er12_1->scale(1.5);

    // Singlet pairs are a bit trickier
    int ij_s = 0;
    for(int i=0; i<nocc_act; i++)
      for(int j=0; j<=i; j++, ij_s++) {
        int ij_ab = i*nocc_act + j;
        int ij_aa = i*(i-1)/2 + j;
        double eab, eaa, e_s;

        eab = epair_ab->get_element(ij_ab);
        if (i != j)
          eaa = epair_aa->get_element(ij_aa);
        else
          eaa = 0.0;
        e_s = (i != j ? 2.0 : 1.0) * eab - 0.5 * eaa;
        epair_0_->set_element(ij_s,e_s);

        eab = emp2_ab->get_element(ij_ab);
        if (i != j)
          eaa = emp2_aa->get_element(ij_aa);
        else
          eaa = 0.0;
        e_s = (i != j ? 2.0 : 1.0) * eab - 0.5 * eaa;
        emp2_0->set_element(ij_s,e_s);

        eab = er12_ab->get_element(ij_ab);
        if (i != j)
          eaa = er12_aa->get_element(ij_aa);
        else
          eaa = 0.0;
        e_s = (i != j ? 2.0 : 1.0) * eab - 0.5 * eaa;
        er12_0->set_element(ij_s,e_s);
      }

    // compute total singlet and triplet energies
    RefSCVector unit_0 = localkit->vector(r12eval_->dim_s());
    RefSCVector unit_1 = localkit->vector(r12eval_->dim_t());
    unit_0->assign(1.0);
    unit_1->assign(1.0);
    emp2tot_0 = emp2_0.dot(unit_0);
    emp2tot_1 = emp2_1.dot(unit_1);
    er12tot_0 = er12_0.dot(unit_0);
    er12tot_1 = er12_1.dot(unit_1);

    ExEnv::out0() << endl << indent << "Singlet MBPT2-R12/" << SA_str << " pair energies:" << endl;
    ExEnv::out0() << indent << scprintf("    i       j        mp2(ij)        r12(ij)      mp2-r12(ij)") << endl;
    ExEnv::out0() << indent << scprintf("  -----   -----   ------------   ------------   ------------") << endl;
    for(int i=0,ij=0;i<nocc_act;i++)
      for(int j=0;j<=i;j++,ij++) {
        ExEnv::out0() << indent << scprintf("  %3d     %3d     %12.9lf   %12.9lf   %12.9lf",i+1,j+1,
                                            emp2_0->get_element(ij),
                                            er12_0->get_element(ij),
                                            epair_0_->get_element(ij)) << endl;
      }
      
    ExEnv::out0() << endl << indent << "Triplet MBPT2-R12/" << SA_str << " pair energies:" << endl;
    ExEnv::out0() << indent << scprintf("    i       j        mp2(ij)        r12(ij)      mp2-r12(ij)") << endl;
    ExEnv::out0() << indent << scprintf("  -----   -----   ------------   ------------   ------------") << endl;
    for(int i=0,ij=0;i<nocc_act;i++)
      for(int j=0;j<i;j++,ij++) {
        ExEnv::out0() << indent << scprintf("  %3d     %3d     %12.9lf   %12.9lf   %12.9lf",i+1,j+1,
                                            emp2_1->get_element(ij),
                                            er12_1->get_element(ij),
                                            epair_1_->get_element(ij)) << endl;
      }

  }


  mp2_corr_energy_ = emp2tot_aa + emp2tot_ab;
  r12_corr_energy_ = er12tot_aa + er12tot_ab;

  delete[] evals;
  tim_exit("mp2-r12a energy");
  
  ///////////////////////////////////////////////////////////////
  // The computation of the MP2 energy is now complete on each
  // node;
  ///////////////////////////////////////////////////////////////

  if (spinadapted_) {
    ExEnv::out0()<<endl<<indent
		 <<scprintf("Singlet MP2 correlation energy [au]:          %17.12lf\n", emp2tot_0);
    ExEnv::out0()<<indent
		 <<scprintf("Triplet MP2 correlation energy [au]:          %17.12lf\n", emp2tot_1);
    ExEnv::out0()<<indent
		 <<scprintf("Singlet (MP2)-R12/%2s correlation energy [au]: %17.12lf\n", SA_str, er12tot_0);
    ExEnv::out0()<<indent
		 <<scprintf("Triplet (MP2)-R12/%2s correlation energy [au]: %17.12lf\n", SA_str, er12tot_1);
    ExEnv::out0()<<indent
		 <<scprintf("Singlet MP2-R12/%2s correlation energy [au]:   %17.12lf\n", SA_str,
			    emp2tot_0 + er12tot_0);
    ExEnv::out0()<<indent
		 <<scprintf("Triplet MP2-R12/%2s correlation energy [au]:   %17.12lf\n", SA_str,
			    emp2tot_1 + er12tot_1);
  }
  
  double etotal = escf + mp2_corr_energy_ + r12_corr_energy_;
  ExEnv::out0()<<endl<<indent
	       <<scprintf("RHF energy [au]:                           %17.12lf\n", escf);
  ExEnv::out0()<<indent
	       <<scprintf("MP2 correlation energy [au]:               %17.12lf\n", mp2_corr_energy_);
  ExEnv::out0()<<indent
	       <<scprintf("(MBPT2)-R12/%2s correlation energy [au]:    %17.12lf\n", SA_str, r12_corr_energy_);
  ExEnv::out0()<<indent
	       <<scprintf("MBPT2-R12/%2s correlation energy [au]:      %17.12lf\n", SA_str, mp2_corr_energy_ + r12_corr_energy_);
  ExEnv::out0()<<indent
	       <<scprintf("MBPT2-R12/%2s energy [au]:                  %17.12lf\n", SA_str, etotal) << endl;

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
