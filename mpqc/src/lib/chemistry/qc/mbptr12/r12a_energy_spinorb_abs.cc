//
// r12a_energy_spinorb_abs.cc
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

#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include <scconfig.h>
#include <util/misc/formio.h>
#include <util/misc/timer.h>
#include <util/group/memory.h>
#include <util/group/message.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/state/state_text.h>
#include <util/state/state_bin.h>
#include <math/scmat/matrix.h>
#include <math/scmat/blocked.h>
#include <math/scmat/repl.h>
#include <math/scmat/local.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/mbpt/util.h>
#include <chemistry/qc/mbpt/bzerofast.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/trans12_grt.h>
#include <chemistry/qc/mbptr12/r12ia.h>
#include <chemistry/qc/mbptr12/r12ia_memgrp.h>
#include <chemistry/qc/mbptr12/r12ia_node0file.h>
#ifdef HAVE_MPIIO
  #include <chemistry/qc/mbptr12/r12ia_mpiiofile.h>
#endif

using namespace std;
using namespace sc;


void
MBPT2_R12::compute_r12a_energy_spinorb_()
{
  double escf = 0.0;
  double emp2 = 0.0;
  double er12a = 0.0;

  tim_enter("mp2-r12a energy");

  if (mp2_done_) {
    escf = hf_energy_;
    emp2 = mp2_corr_energy_;
  }
  else {
    escf = ref_energy();
    emp2 = MBPT2::corr_energy();
  }

  // Convert V's into local SCMatrices on node 0
  // Then compute B (local SCMatrices as well)
  Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();

  int nocc = 0;
  for (int i=0; i<oso_dimension()->n(); i++) {
    if (reference_->occupation(i) == 2.0) nocc++;
  }
  int nocc_act = nocc - nfzc;
  int naa = nocc_act*(nocc_act-1)/2;
  int nab = nocc_act*nocc_act;
  RefSCDimension npair_ab = new SCDimension(nocc_act*nocc_act);
  RefSCDimension npair_aa = new SCDimension((nocc_act*(nocc_act-1))/2);
  int me = msg_->me();
  if (me == 0) {
    Vaa = local_matrix_kit->matrix(npair_aa,npair_aa);
    Xaa = local_matrix_kit->matrix(npair_aa,npair_aa);
    Baa = local_matrix_kit->matrix(npair_aa,npair_aa);
    Vab = local_matrix_kit->matrix(npair_ab,npair_ab);
    Xab = local_matrix_kit->matrix(npair_ab,npair_ab);
    Bab = local_matrix_kit->matrix(npair_ab,npair_ab);
    Vaa->unit();
    Vab->unit();
    Baa->unit();
    Bab->unit();
    Xaa->assign(0.0);
    Xab->assign(0.0);

    if (debug_)
      ExEnv::out0() << indent << "Allocated V, X, and B" << endl;

    compute_r12a_intermed_spinorb_SBS_();
    compute_r12a_intermed_spinorb_ABS_();

    if (debug_) {
      Vaa->print("Alpha-alpha V matrix");
      Baa->print("Alpha-alpha B matrix");
      Vab->print("Alpha-beta V matrix");
      Bab->print("Alpha-beta B matrix");
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

    //
    // Evaluate pair energies
    //
    epair_aa = local_matrix_kit->vector(npair_aa);
    epair_ab = local_matrix_kit->vector(npair_ab);

    // For some reason invert_this doesn't work here
    Baa->gen_invert_this();
    Bab->gen_invert_this();
    if (debug_ > 1) {
      Baa->print("Inverse alpha-alpha B matrix");
      Bab->print("Inverse alpha-beta B matrix");
    }
  
    double er12a_aa = 0.0;
    double er12a_ab = 0.0;
    
    ExEnv::out0() << endl << indent << "Alpha-alpha MBPT2-R12/A pair energies:" << endl;
    ExEnv::out0() << indent << scprintf("    i       j         e(ij)") << endl;
    ExEnv::out0() << indent << scprintf("  -----   -----   ------------") << endl;
    int ij, kl;
    for(int i=0,ij=0;i<nocc_act;i++)
      for(int j=0;j<i;j++,ij++) {
	double eij = 0.0;
	for(kl=0;kl<naa;kl++)
	  for(int mn=0;mn<naa;mn++) {
	    eij -= 2.0*Vaa->get_element(kl,ij)*Vaa->get_element(mn,ij)*Baa->get_element(kl,mn);
	  }
	epair_aa->set_element(ij,eij);
	er12a_aa += eij;
	ExEnv::out0() << indent << scprintf("  %3d     %3d     %12.9lf",i+1,j+1,eij) << endl;
      }
    
    ExEnv::out0() << endl << indent << "Alpha-beta MBPT2-R12/A pair energies:" << endl;
    ExEnv::out0() << indent << scprintf("    i       j         e(ij)") << endl;
    ExEnv::out0() << indent << scprintf("  -----   -----   ------------") << endl;
    for(int i=0,ij=0;i<nocc_act;i++)
      for(int j=0;j<nocc_act;j++,ij++) {
	double eij = 0.0;
	for(kl=0;kl<nab;kl++)
	  for(int mn=0;mn<nab;mn++) {
	    eij -= 1.0*Vab->get_element(kl,ij)*Vab->get_element(mn,ij)*Bab->get_element(kl,mn);
	  }
	epair_ab->set_element(ij,eij);
	er12a_ab += eij;
	ExEnv::out0() << indent << scprintf("  %3d     %3d     %12.9lf",i+1,j+1,eij) << endl;
      }

    er12a = er12a_aa + er12a_ab;
    r12_corr_energy_ = er12a;
  }
  tim_exit("mp2-r12a energy");

  ///////////////////////////////////////////////////////////////
  // The computation of the MP2 energy is now complete on each
  // node;
  ///////////////////////////////////////////////////////////////

  double emp2r12a = escf + emp2 + er12a;
  ExEnv::out0()<<endl<<indent
	       <<scprintf("RHF energy [au]:                           %17.12lf\n", escf);
  ExEnv::out0()<<indent
	       <<scprintf("MP2 correlation energy [au]:               %17.12lf\n", emp2);
  ExEnv::out0()<<indent
	       <<scprintf("(MBPT2)-R12/A correlation energy [au]:     %17.12lf\n", er12a);
  ExEnv::out0()<<indent
	       <<scprintf("MBPT2-R12/A energy [au]:                   %17.12lf\n", emp2r12a);
  ExEnv::out0().flush();

  set_energy(emp2r12a);
  set_actual_value_accuracy(reference_->actual_value_accuracy()
                            *ref_to_mp2_acc);

  return;
}


////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
