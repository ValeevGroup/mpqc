//
// dualbasis_mp2.cc
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

#include <stdexcept>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include <scconfig.h>
#include <util/misc/formio.h>
#include <util/misc/regtime.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/state/state_text.h>
#include <util/state/state_bin.h>
#include <math/scmat/matrix.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/mbpt/bzerofast.h>
#include <chemistry/qc/mbptr12/r12ia.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>
#include <chemistry/qc/mbptr12/pairiter.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>

using namespace std;
using namespace sc;

void
R12IntEval::compute_dualEmp2_()
{
  if (evaluated_)
    return;
  Ref<MessageGrp> msg = r12info()->msg();
  Ref<MemoryGrp> mem = r12info()->mem();
  Ref<ThreadGrp> thr = r12info()->thr();
  const int num_te_types = 1;
  enum te_types {eri=0};

  Timer tim("dual-basis MP2 energy");
  ExEnv::out0() << endl << indent
	       << "Entered dual-basis MP2 energy evaluator" << endl;
  ExEnv::out0() << incindent;

  int me = msg->me();
  int nproc = msg->n();
  
  // Do the AO->MO transform
  form_canonvir_space_();
  Ref<MOIntsTransformFactory> tfactory = r12info_->tfactory();
  tfactory->set_spaces(r12info_->act_occ_space(),canonvir_space_,
                       r12info_->act_occ_space(),canonvir_space_);
  Ref<TwoBodyMOIntsTransform> ipjq_tform = tfactory->twobody_transform_13("(ix|jy)");
  ipjq_tform->set_num_te_types(num_te_types);
  ipjq_tform->compute();
  Ref<R12IntsAcc> ijpq_acc = ipjq_tform->ints_acc();
  if (num_te_types != ijpq_acc->num_te_types())
    throw std::runtime_error("R12IntEval::compute_dualEmp2_() -- number of MO integral types is wrong");

  int nocc_act = r12info()->nocc_act();
  int ncanonvir = canonvir_space_->rank();

  ExEnv::out0() << indent << "Begin computation of energies" << endl;
  SpatialMOPairIter_eq kl_iter(r12info_->act_occ_space());
  int naa = kl_iter.nij_aa();          // Number of alpha-alpha pairs (i > j)
  int nab = kl_iter.nij_ab();          // Number of alpha-beta pairs
  if (debug_) {
    ExEnv::out0() << indent << "naa = " << naa << endl;
    ExEnv::out0() << indent << "nab = " << nab << endl;
  }

  // Compute the number of tasks that have full access to the integrals
  // and split the work among them
  vector<int> proc_with_ints;
  int nproc_with_ints = tasks_with_ints_(ijpq_acc,proc_with_ints);

  RefDiagSCMatrix act_occ_evals = r12info_->act_occ_space()->evals();
  RefDiagSCMatrix canonvir_evals = canonvir_space_->evals();

  if (ijpq_acc->has_access(me)) {

    for(kl_iter.start();int(kl_iter);kl_iter.next()) {

      const int kl = kl_iter.ij();
      // Figure out if this task will handle this kl
      int kl_proc = kl%nproc_with_ints;
      if (kl_proc != proc_with_ints[me])
        continue;
      const int k = kl_iter.i();
      const int l = kl_iter.j();
      const int kl_aa = kl_iter.ij_aa();
      const int kl_ab = kl_iter.ij_ab();
      const int lk_ab = kl_iter.ij_ba();

      if (debug_)
        ExEnv::outn() << indent << "task " << me << ": working on (k,l) = " << k << "," << l << " " << endl;

      // Get (|1/r12|) integrals
      tim.enter("MO ints retrieve");
      double *klxy_buf_eri = ijpq_acc->retrieve_pair_block(k,l,R12IntsAcc::eri);
      tim.exit("MO ints retrieve");

      if (debug_)
        ExEnv::outn() << indent << "task " << me << ": obtained kl blocks" << endl;

      // Compute MP2 energies
      double emp2_aa = 0.0;
      double emp2_ab = 0.0;
      for(int a=0; a<ncanonvir; a++) {
        for(int b=0; b<ncanonvir; b++) {
          const int ab_offset = a*ncanonvir+b;
          const int ba_offset = b*ncanonvir+a;
          const double oo_delta_ijab = 1.0/(act_occ_evals(k)+act_occ_evals(l)-canonvir_evals(a)-canonvir_evals(b));
          const double eri_kalb = klxy_buf_eri[ab_offset];
          const double eri_kbla = klxy_buf_eri[ba_offset];
          emp2_ab += 0.5*(eri_kalb * eri_kalb + eri_kbla * eri_kbla) * oo_delta_ijab;
          if (kl_aa != -1) {
            emp2_aa += (eri_kalb - eri_kbla) * (eri_kalb - eri_kbla) * oo_delta_ijab;
          }
        }
      }
      emp2pair_ab_.set_element(kl_ab,emp2_ab);
      if (kl_ab != lk_ab)
        emp2pair_ab_.set_element(lk_ab,emp2_ab);
      if (kl_aa != -1)
        emp2pair_aa_.set_element(kl_aa,emp2_aa);

      ijpq_acc->release_pair_block(k,l,R12IntsAcc::eri);
    }
  }

  // Tasks that don't do any work here still need to create these timers
  tim.enter("MO ints retrieve");
  tim.exit("MO ints retrieve");

  ExEnv::out0() << indent << "End of computation of energies" << endl;
  ijpq_acc->deactivate();
  
  globally_sum_intermeds_();

  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Exited dual-basis MP2 energy evaluator" << endl;

  tim.exit("dual-basis MP2 energy");
  checkpoint_();
  
  return;
}

void
R12IntEval::compute_dualEmp1_()
{
  if (evaluated_)
    return;
  Ref<MessageGrp> msg = r12info()->msg();
  Ref<MemoryGrp> mem = r12info()->mem();
  Ref<ThreadGrp> thr = r12info()->thr();
  const int num_te_types = 1;
  enum te_types {eri=0};

  Timer tim("dual-basis MP1 energy");
  ExEnv::out0() << endl << indent
	       << "Entered dual-basis MP1 energy evaluator" << endl;
  ExEnv::out0() << incindent;

  int me = msg->me();
  int nproc = msg->n();
  
  // Compute act.occ./aux.virt. Fock matrix
  form_canonvir_space_();
  Ref<MOIndexSpace> occ_space = r12info_->occ_space();
  RefSCMatrix F_aocc_canonvir = fock_(occ_space,occ_space,canonvir_space_);

  int nocc = r12info()->nocc();
  int ncanonvir = canonvir_space_->rank();
  RefDiagSCMatrix occ_evals = r12info_->occ_space()->evals();
  RefDiagSCMatrix canonvir_evals = canonvir_space_->evals();

  double emp1 = 0.0;
  for(int i=0; i<nocc; i++) {
    for(int a=0; a<ncanonvir; a++) {
      const double Fia = F_aocc_canonvir.get_element(i,a);
      emp1 += Fia*Fia/(-occ_evals(i)+canonvir_evals(a));
    }
  }
  ExEnv::out0() << indent << "MP1 energy correction to HF energy [au] :   "
                << 2.0*emp1 << endl;
  ExEnv::out0() << indent << "HF energy estimated in new basis [au]   :   "
                << r12info_->ref()->energy() - 2.0*emp1 << endl;

  ExEnv::out0() << decindent;
  ExEnv::out0() << endl << "Exited dual-basis MP1 energy evaluator" << endl;

  tim.exit("dual-basis MP1 energy");
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
