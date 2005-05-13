//
// intermeds_g12.cc
//
// Copyright (C) 2005 Edward Valeev
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
#include <util/misc/timer.h>
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
R12IntEval::init_intermeds_g12_()
{
  Ref<MessageGrp> msg = r12info()->msg();
  Ref<MemoryGrp> mem = r12info()->mem();
  Ref<ThreadGrp> thr = r12info()->thr();

  tim_enter("\"diagonal\" part of G12 intermediates");
  ExEnv::out0() << endl << indent
	       << "Entered G12 diagonal intermediates evaluator" << endl;
  ExEnv::out0() << incindent;

  int me = msg->me();
  int nproc = msg->n();

  //
  // Do the AO->MO transform:
  // 1) get (im|jn) integrals of g12/r12 and [g12,[t1,g12]] operators to compute diagonal parts of V and B
  //
  Ref<MOIntsTransformFactory> tfactory = r12info_->tfactory();
  Ref<TwoBodyMOIntsTransform> imjn_tform = tform_map_["(im|jn)"];
  Ref<R12IntsAcc> ijmn_acc = imjn_tform->ints_acc();
  if (!ijmn_acc->is_committed()) {
    imjn_tform->compute(corrparam_);
    }

  //
  // 2) get (im|jn) integrals of g12*g12 operator (in reality use g12 integrals with the exponent multiplied by 2)
  //    these integrals used to compute X
  //    NOTE: use occ_space instead of act_occ_space so that one block of code will handle all 3 integrals
  tfactory->set_spaces(r12info_->act_occ_space(),r12info->occ_space(),
                       r12info_->act_occ_space(),r12info->occ_space());
  Ref<TwoBodyMOIntsTransform> im2jn_tform = tfactory->twobody_transform_13("(im|2|jn)",corrfactor_->callback());
  im2jn_tform->compute(2.0*corrparam_);
  Ref<R12IntsAcc> ij2mn_acc = im2jn_tform->ints_acc();

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
  int nproc_with_ints = tasks_with_ints_(ijmn_acc,proc_with_ints);

  if (ijmn_acc->has_access(me)) {

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
      tim_enter("MO ints retrieve");
      double *klxy_buf_eri = ijmn_acc->retrieve_pair_block(k,l,corrfactor_->tbint_type_eri());
      tim_exit("MO ints retrieve");

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

      ijmn_acc->release_pair_block(k,l,corrfactor_->tbint_type_eri());
    }
  }

  // Tasks that don't do any work here still need to create these timers
  tim_enter("MO ints retrieve");
  tim_exit("MO ints retrieve");

  ExEnv::out0() << indent << "End of computation of energies" << endl;
  ijmn_acc->deactivate();
  
  globally_sum_intermeds_();

  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Exited dual-basis MP2 energy evaluator" << endl;

  tim_exit("dual-basis MP2 energy");
  checkpoint_();
  
  return;
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
