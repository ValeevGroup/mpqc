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
  // 1) get (im|jn) integrals of g12/r12 operator to compute diagonal parts of V
  //
  Ref<MOIntsTransformFactory> tfactory = r12info_->tfactory();
  Ref<TwoBodyMOIntsTransform> imjn_tform = tform_map_["(im|jn)"];
  Ref<R12IntsAcc> ijmn_acc = imjn_tform->ints_acc();
  if (!ijmn_acc->is_committed()) {
    imjn_tform->compute(corrparam_);
  }

  //
  // 2) get (im|jn) integrals of [g12,[t1,g12]] and g12*g12 operator (use integrals with the exponent multiplied by 2, and additionally [g12,[t1,g12]] integral needs to be scaled by 0.25 to take into account that real exponent is half what the integral library thinks)
  //    these integrals used to compute X and B
  //    NOTE: use occ_space instead of act_occ_space so that one block of code will handle all 3 integrals
  tfactory->set_spaces(r12info_->act_occ_space(),r12info_->occ_space(),
                       r12info_->act_occ_space(),r12info_->occ_space());
  Ref<TwoBodyMOIntsTransform> im2jn_tform = tfactory->twobody_transform_13("(im|2|jn)",corrfactor_->callback());
  im2jn_tform->set_num_te_types(corrfactor_->num_tbint_types());
  im2jn_tform->compute(2.0*corrparam_);
  Ref<R12IntsAcc> ij2mn_acc = im2jn_tform->ints_acc();

  int nfzc = r12info()->nfzc();
  int nocc_act = r12info()->ndocc_act();
  int nocc = r12info()->ndocc();

  ExEnv::out0() << indent << "Begin computation of intermediates" << endl;
  SpatialMOPairIter_eq ij_iter(r12info_->act_occ_space());
  int naa = ij_iter.nij_aa();          // Number of alpha-alpha pairs (i > j)
  int nab = ij_iter.nij_ab();          // Number of alpha-beta pairs
  if (debug_) {
    ExEnv::out0() << indent << "naa = " << naa << endl;
    ExEnv::out0() << indent << "nab = " << nab << endl;
  }

  // Compute the number of tasks that have full access to the integrals
  // and split the work among them
  vector<int> proc_with_ints;
  int nproc_with_ints = tasks_with_ints_(ijmn_acc,proc_with_ints);

  if (ijmn_acc->has_access(me)) {

    for(ij_iter.start();int(ij_iter);ij_iter.next()) {

      const int ij = ij_iter.ij();
      // Figure out if this task will handle this ij
      int ij_proc = ij%nproc_with_ints;
      if (ij_proc != proc_with_ints[me])
        continue;
      const int i = ij_iter.i();
      const int j = ij_iter.j();
      const int ij_aa = ij_iter.ij_aa();
      const int ij_ab = ij_iter.ij_ab();
      const int ji_ab = ij_iter.ij_ba();

      if (debug_)
        ExEnv::outn() << indent << "task " << me << ": working on (i,j) = " << i << "," << j << " " << endl;

      // Get the integrals
      tim_enter("MO ints retrieve");
      double *ijxy_buf_f12eri   = ijmn_acc->retrieve_pair_block(i,j,corrfactor_->tbint_type_f12eri());
      double *ijxy_buf_f12t1f12 = ij2mn_acc->retrieve_pair_block(i,j,corrfactor_->tbint_type_f12t1f12());
      double *ijxy_buf_f12f12   = ij2mn_acc->retrieve_pair_block(i,j,corrfactor_->tbint_type_f12f12());
      tim_exit("MO ints retrieve");

      if (debug_)
        ExEnv::outn() << indent << "task " << me << ": obtained ij blocks" << endl;

      SpatialMOPairIter_eq kl_iter(r12info_->act_occ_space());
      for(kl_iter.start();int(kl_iter);kl_iter.next()) {
        const int k = kl_iter.i();
        const int l = kl_iter.j();
        const int kl_aa = kl_iter.ij_aa();
        const int kl_ab = kl_iter.ij_ab();
        const int lk_ab = kl_iter.ij_ba();
        
        const int kk = k + nfzc;
        const int ll = l + nfzc;
        const int kkll = kk*nocc+ll;
        const int llkk = ll*nocc+kk;
        
        const double V_ijkl_ab = ijxy_buf_f12eri[kkll];
        const double V_ijlk_ab = ijxy_buf_f12eri[llkk];
        const double V_jikl_ab = V_ijlk_ab;
        const double V_jilk_ab = V_ijkl_ab;
        Vab_.set_element(ij_ab,kl_ab,V_ijkl_ab);
        Vab_.set_element(ji_ab,kl_ab,V_jikl_ab);
        Vab_.set_element(ij_ab,lk_ab,V_ijlk_ab);
        Vab_.set_element(ji_ab,lk_ab,V_jilk_ab);
        if (ij_aa != -1 && kl_aa != -1)
          Vaa_.set_element(ij_aa,kl_aa,V_ijkl_ab-V_ijlk_ab);
      
        const double X_ijkl_ab = ijxy_buf_f12f12[kkll];
        const double X_ijlk_ab = ijxy_buf_f12f12[llkk];
        const double X_jikl_ab = X_ijlk_ab;
        const double X_jilk_ab = X_ijkl_ab;
        Xab_.set_element(ij_ab,kl_ab,X_ijkl_ab);
        Xab_.set_element(ji_ab,kl_ab,X_jikl_ab);
        Xab_.set_element(ij_ab,lk_ab,X_ijlk_ab);
        Xab_.set_element(ji_ab,lk_ab,X_jilk_ab);
        if (ij_aa != -1 && kl_aa != -1)
          Xaa_.set_element(ij_aa,kl_aa,X_ijkl_ab-X_ijlk_ab);

        const double B_ijkl_ab = 0.25 * ijxy_buf_f12t1f12[kkll];
        const double B_ijlk_ab = 0.25 * ijxy_buf_f12t1f12[llkk];
        const double B_jikl_ab = B_ijlk_ab;
        const double B_jilk_ab = B_ijkl_ab;
        Bab_.set_element(ij_ab,kl_ab,B_ijkl_ab);
        Bab_.set_element(ji_ab,kl_ab,B_jikl_ab);
        Bab_.set_element(ij_ab,lk_ab,B_ijlk_ab);
        Bab_.set_element(ji_ab,lk_ab,B_jilk_ab);
        if (ij_aa != -1 && kl_aa != -1)
          Baa_.set_element(ij_aa,kl_aa,B_ijkl_ab-B_ijlk_ab);
        
      }
      
      ij2mn_acc->release_pair_block(i,j,corrfactor_->tbint_type_f12f12());
      ij2mn_acc->release_pair_block(i,j,corrfactor_->tbint_type_f12t1f12());
      ijmn_acc->release_pair_block(i,j,corrfactor_->tbint_type_f12eri());
    }
  }

  // Tasks that don't do any work here still need to create these timers
  tim_enter("MO ints retrieve");
  tim_exit("MO ints retrieve");

  ExEnv::out0() << indent << "End of computation of intermediates" << endl;
  ijmn_acc->deactivate();
  ij2mn_acc->deactivate();
  
  globally_sum_intermeds_();

  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Exited G12 diagonal intermediates evaluator" << endl;

  tim_exit("\"diagonal\" part of G12 intermediates");
  checkpoint_();
  
  return;
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
