//
// compute_vxb_a_symm.cc
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
#include <sstream>
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

#define PRINT_R12_INTERMED 0

void
R12IntEval::contrib_to_VXB_a_symm_(const std::string& tform_name)
{
  if (evaluated_)
    return;
  LinearR12::ABSMethod abs_method = r12info_->abs_method();
  Ref<MessageGrp> msg = r12info_->msg();
  Ref<MemoryGrp> mem = r12info_->mem();
  Ref<ThreadGrp> thr = r12info_->thr();

  tim_enter("mp2-r12a intermeds (symmetric term)");

  int me = msg->me();
  int nproc = msg->n();
  
  // Do the AO->MO transform
  Ref<TwoBodyMOIntsTransform> ipjq_tform = get_tform_(tform_name);
  if (ipjq_tform->space1() != ipjq_tform->space3())
    throw std::runtime_error("R12IntEval::contrib_to_VXB_a_symm_() -- wrong type of transform is provided (space1 != space3)");
  Ref<MOIndexSpace> mospace = ipjq_tform->space2();

  // Carry out the AO->MO transform
  Ref<R12IntsAcc> ijpq_acc = ipjq_tform->ints_acc();
  if (ijpq_acc.null() || !ijpq_acc->is_committed())
    ipjq_tform->compute(intparams_);
  if (!ijpq_acc->is_active())
    ijpq_acc->activate();

  ExEnv::out0() << endl << indent
                << "Entered \"" << mospace->name()
                << "\" A (GEBC) intermediates evaluator" << endl;
  ExEnv::out0() << incindent;

  const Ref<MOIndexSpace>& act_occ_space = ipjq_tform->space1();
  const int nocc_act = act_occ_space->rank();
  const int noso = mospace->rank();

  /*--------------------------------
    Compute MP2-R12/A intermediates
    and collect on node0
   --------------------------------*/
  ExEnv::out0() << indent << "Begin computation of intermediates" << endl;
  tim_enter("intermediates");
  SpatialMOPairIter_eq ij_iter(act_occ_space);
  SpatialMOPairIter_eq kl_iter(act_occ_space);
  int naa = ij_iter.nij_aa();          // Number of alpha-alpha pairs (i > j)
  int nab = ij_iter.nij_ab();          // Number of alpha-beta pairs
  if (debug_) {
    ExEnv::out0() << indent << "naa = " << naa << endl;
    ExEnv::out0() << indent << "nab = " << nab << endl;
  }

  // Compute intermediates
  if (debug_)
    ExEnv::out0() << indent << "Ready to compute MP2-R12/A (GEBC) intermediates" << endl;

  // Compute the number of tasks that have full access to the integrals
  // and split the work among them
  vector<int> proc_with_ints;
  const int nproc_with_ints = tasks_with_ints_(ijpq_acc,proc_with_ints);

  
  //////////////////////////////////////////////////////////////
  //
  // Evaluation of the intermediates proceeds as follows:
  //
  //    loop over batches of kl, k >= l,  0<=k,l<nocc_act
  //      load (kl|xy), (kl| [T1,r12] |xy), and (lk| [T1,r12] |xy)
  //           (aka kl-sets) into memory
  //
  //      loop over batches of ij, i>=j, 0<=i,j<nocc_act
  //        load (ij|r12|xy) into memory
  //           (aka ij-sets) into memory
  //        compute V[ij][kl] and T[ij][kl] for all ij and kl in
  //                the "direct product" batch
  //      end ij loop
  //    end kl loop
  //
  /////////////////////////////////////////////////////////////////////////////////

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

      // Get (|1/r12|), (|r12|), and (|[r12,T1]|) integrals
      tim_enter("MO ints retrieve");
      double *klxy_buf_eri = ijpq_acc->retrieve_pair_block(k,l,corrfactor_->tbint_type_eri());
      double *klxy_buf_r12 = ijpq_acc->retrieve_pair_block(k,l,corrfactor_->tbint_type_f12());
      double *klxy_buf_r12t1 = ijpq_acc->retrieve_pair_block(k,l,corrfactor_->tbint_type_t1f12());
      double *lkxy_buf_r12t1 = ijpq_acc->retrieve_pair_block(l,k,corrfactor_->tbint_type_t1f12());
      tim_exit("MO ints retrieve");

      if (debug_)
        ExEnv::outn() << indent << "task " << me << ": obtained kl blocks" << endl;

      // to avoid every task hitting same ij at the same time, stagger ij-accesses, i.e. each kl task will start with ij=kl+1
      for(ij_iter.start(kl+1);int(ij_iter);ij_iter.next()) {

        const int i = ij_iter.i();
        const int j = ij_iter.j();
        const int ij_aa = ij_iter.ij_aa();
        const int ij_ab = ij_iter.ij_ab();
        const int ji_ab = ij_iter.ij_ba();

        if (debug_)
          ExEnv::outn() << indent << "task " << me << ": (k,l) = " << k << "," << l << ": (i,j) = " << i << "," << j << endl;

        tim_enter("MO ints retrieve");
        double *ijxy_buf_r12 = ijpq_acc->retrieve_pair_block(i,j,corrfactor_->tbint_type_f12());
        tim_exit("MO ints retrieve");

        if (debug_)
          ExEnv::outn() << indent << "task " << me << ": obtained ij blocks" << endl;


        tim_enter("MO ints contraction");
        double Vaa_ijkl, Vab_ijkl, Vab_jikl, Vab_ijlk, Vab_jilk;
        double Xaa_ijkl, Xab_ijkl, Xab_jikl, Xab_ijlk, Xab_jilk;
        double Taa_ijkl, Tab_ijkl, Tab_jikl, Tab_ijlk, Tab_jilk;
        Vaa_ijkl = Vab_ijkl = Vab_jikl = Vab_ijlk = Vab_jilk = 0.0;
        Xaa_ijkl = Xab_ijkl = Xab_jikl = Xab_ijlk = Xab_jilk = 0.0;
        Taa_ijkl = Tab_ijkl = Tab_jikl = Tab_ijlk = Tab_jilk = 0.0;

        const double pfac_xy = 0.5;
        for(int y=0;y<noso;y++) {
          for(int x=0;x<noso;x++) {
            int yx_offset = y*noso+x;
            int xy_offset = x*noso+y;
            double ij_r12_xy = ijxy_buf_r12[xy_offset];
            double ij_r12_yx = ijxy_buf_r12[yx_offset];
            double kl_eri_xy = klxy_buf_eri[xy_offset];
            double kl_eri_yx = klxy_buf_eri[yx_offset];
            Vab_ijkl -= pfac_xy * (ij_r12_xy * kl_eri_xy + ij_r12_yx * kl_eri_yx);
            if (ij_ab != ji_ab)
              Vab_jikl -= pfac_xy * (ij_r12_yx * kl_eri_xy + ij_r12_xy * kl_eri_yx);
            if (kl_ab != lk_ab)
              Vab_ijlk -= pfac_xy * (ij_r12_xy * kl_eri_yx + ij_r12_yx * kl_eri_xy);
            if (ij_ab != ji_ab && kl_ab != lk_ab) {
              Vab_jilk -= pfac_xy * (ij_r12_yx * kl_eri_yx + ij_r12_xy * kl_eri_xy);
            }
            if (ij_aa != -1 && kl_aa != -1) {
              Vaa_ijkl -= pfac_xy * (ij_r12_xy - ij_r12_yx)*(kl_eri_xy - kl_eri_yx);
            }
            double kl_r12_xy = klxy_buf_r12[xy_offset];
            double kl_r12_yx = klxy_buf_r12[yx_offset];
            Xab_ijkl -= pfac_xy * (ij_r12_xy * kl_r12_xy + ij_r12_yx * kl_r12_yx);
            if (ij_ab != ji_ab)
              Xab_jikl -= pfac_xy * (ij_r12_yx * kl_r12_xy + ij_r12_xy * kl_r12_yx);
            if (kl_ab != lk_ab)
              Xab_ijlk -= pfac_xy * (ij_r12_xy * kl_r12_yx + ij_r12_yx * kl_r12_xy);
            if (ij_ab != ji_ab && kl_ab != lk_ab) {
              Xab_jilk -= pfac_xy * (ij_r12_yx * kl_r12_yx + ij_r12_xy * kl_r12_xy);
            }
            if (ij_aa != -1 && kl_aa != -1) {
              Xaa_ijkl -= pfac_xy * (ij_r12_xy - ij_r12_yx)*(kl_r12_xy - kl_r12_yx);
            }
            double kl_r12t1_xy = klxy_buf_r12t1[xy_offset];
            double kl_r12t1_yx = klxy_buf_r12t1[yx_offset];
            double lk_r12t1_xy = lkxy_buf_r12t1[xy_offset];
            double lk_r12t1_yx = lkxy_buf_r12t1[yx_offset];
            double kl_Tr12_xy = -kl_r12t1_xy-lk_r12t1_yx;
            double kl_Tr12_yx = -kl_r12t1_yx-lk_r12t1_xy;
            Tab_ijkl += pfac_xy * (ij_r12_xy * kl_Tr12_xy + ij_r12_yx * kl_Tr12_yx);
            if (ij_ab != ji_ab)
              Tab_jikl += pfac_xy * (ij_r12_yx * kl_Tr12_xy + ij_r12_xy * kl_Tr12_yx);
            if (kl_ab != lk_ab)
              Tab_ijlk += pfac_xy * (ij_r12_xy * kl_Tr12_yx + ij_r12_yx * kl_Tr12_xy);
            if (ij_ab != ji_ab && kl_ab != lk_ab) {
              Tab_jilk += pfac_xy * (ij_r12_yx * kl_Tr12_yx + ij_r12_xy * kl_Tr12_xy);
            }
            if (ij_aa != -1 && kl_aa != -1) {
              Taa_ijkl += pfac_xy * (ij_r12_xy - ij_r12_yx)*(kl_Tr12_xy - kl_Tr12_yx);
            }
          }
        }
        Vab_.accumulate_element(ij_ab,kl_ab,Vab_ijkl);
        if (ij_ab != ji_ab)
          Vab_.accumulate_element(ji_ab,kl_ab,Vab_jikl);
        if (kl_ab != lk_ab)
          Vab_.accumulate_element(ij_ab,lk_ab,Vab_ijlk);
        if (ij_ab != ji_ab && kl_ab != lk_ab)
          Vab_.accumulate_element(ji_ab,lk_ab,Vab_jilk);
        if (ij_aa != -1 && kl_aa != -1)
          Vaa_.accumulate_element(ij_aa,kl_aa,Vaa_ijkl);
        Xab_.accumulate_element(ij_ab,kl_ab,Xab_ijkl);
        if (ij_ab != ji_ab)
          Xab_.accumulate_element(ji_ab,kl_ab,Xab_jikl);
        if (kl_ab != lk_ab)
          Xab_.accumulate_element(ij_ab,lk_ab,Xab_ijlk);
        if (ij_ab != ji_ab && kl_ab != lk_ab)
          Xab_.accumulate_element(ji_ab,lk_ab,Xab_jilk);
        if (ij_aa != -1 && kl_aa != -1)
          Xaa_.accumulate_element(ij_aa,kl_aa,Xaa_ijkl);
        Bab_.accumulate_element(ij_ab,kl_ab,Tab_ijkl);
        if (ij_ab != ji_ab)
          Bab_.accumulate_element(ji_ab,kl_ab,Tab_jikl);
        if (kl_ab != lk_ab)
          Bab_.accumulate_element(ij_ab,lk_ab,Tab_ijlk);
        if (ij_ab != ji_ab && kl_ab != lk_ab)
          Bab_.accumulate_element(ji_ab,lk_ab,Tab_jilk);
        if (ij_aa != -1 && kl_aa != -1)
          Baa_.accumulate_element(ij_aa,kl_aa,Taa_ijkl);
        tim_exit("MO ints contraction");

#if PRINT_R12_INTERMED
	    if (ij_ab != ji_ab && kl_ab != lk_ab)
	      printf("Vaa[%d][%d] = %lf\n",ij_aa,kl_aa,Vaa_ij[kl_aa]);
            printf("Vab[%d][%d] = %lf\n",ij_ab,kl_ab,Vab_ij[kl_ab]);
	    if (ij_ab != ji_ab)
	      printf("Vab[%d][%d] = %lf\n",ji_ab,kl_ab,Vab_ji[kl_ab]);
	    if (kl_ab != lk_ab)
	      printf("Vab[%d][%d] = %lf\n",ij_ab,lk_ab,Vab_ij[lk_ab]);
	    if (ij_ab != ji_ab && kl_ab != lk_ab)
	      printf("Vab[%d][%d] = %lf\n",ji_ab,lk_ab,Vab_ji[lk_ab]);
	    if (ij_ab != ji_ab && kl_ab != lk_ab)
	      printf("Xaa[%d][%d] = %lf\n",ij_aa,kl_aa,Xaa_ij[kl_aa]);
            printf("Xab[%d][%d] = %lf\n",ij_ab,kl_ab,Xab_ij[kl_ab]);
	    if (ij_ab != ji_ab)
	      printf("Xab[%d][%d] = %lf\n",ji_ab,kl_ab,Xab_ji[kl_ab]);
	    if (kl_ab != lk_ab)
	      printf("Xab[%d][%d] = %lf\n",ij_ab,lk_ab,Xab_ij[lk_ab]);
	    if (ij_ab != ji_ab && kl_ab != lk_ab)
	      printf("Xab[%d][%d] = %lf\n",ji_ab,lk_ab,Xab_ji[lk_ab]);
	    if (ij_ab != ji_ab && kl_ab != lk_ab)
	      printf("Taa[%d][%d] = %lf\n",ij_aa,kl_aa,Taa_ij[kl_aa]);
            printf("Tab[%d][%d] = %lf\n",ij_ab,kl_ab,Tab_ij[kl_ab]);
	    if (ij_ab != ji_ab)
	      printf("Tab[%d][%d] = %lf\n",ji_ab,kl_ab,Tab_ji[kl_ab]);
	    if (kl_ab != lk_ab)
	      printf("Tab[%d][%d] = %lf\n",ij_ab,lk_ab,Tab_ij[lk_ab]);
	    if (ij_ab != ji_ab && kl_ab != lk_ab)
	      printf("Tab[%d][%d] = %lf\n",ji_ab,lk_ab,Tab_ji[lk_ab]);
#endif
        ijpq_acc->release_pair_block(i,j,corrfactor_->tbint_type_f12());
      }
      ijpq_acc->release_pair_block(k,l,corrfactor_->tbint_type_eri());
      ijpq_acc->release_pair_block(k,l,corrfactor_->tbint_type_f12());
      ijpq_acc->release_pair_block(k,l,corrfactor_->tbint_type_t1f12());
      ijpq_acc->release_pair_block(l,k,corrfactor_->tbint_type_t1f12());
    }
  }
  // Tasks that don't do any work here still need to create these timers
  tim_enter("MO ints retrieve");
  tim_exit("MO ints retrieve");
  tim_enter("MO ints contraction");
  tim_exit("MO ints contraction");

  tim_exit("intermediates");
  ExEnv::out0() << indent << "End of computation of intermediates" << endl;
  ijpq_acc->deactivate();
  
  // Symmetrize B intermediate
  for(int ij=0;ij<naa;ij++)
    for(int kl=0;kl<=ij;kl++) {
      double belem = 0.5*(Baa_->get_element(ij,kl) + Baa_->get_element(kl,ij));
      Baa_->set_element(ij,kl,belem);
      Baa_->set_element(kl,ij,belem);
    }

  for(int ij=0;ij<nab;ij++)
    for(int kl=0;kl<=ij;kl++) {
      double belem = 0.5*(Bab_->get_element(ij,kl) + Bab_->get_element(kl,ij));
      Bab_->set_element(ij,kl,belem);
      Bab_->set_element(kl,ij,belem);
    }

  globally_sum_intermeds_();
  
  ExEnv::out0() << decindent;
  ExEnv::out0() << indent
                << "Exited \"" << mospace->name()
                << "\" A (GEBC) intermediates evaluator" << endl;

  tim_exit("mp2-r12a intermeds (symmetric term)");
  checkpoint_();
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
