//
// compute_a_gebc.cc
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

#define PRINT4Q_MP2 0
#define PRINT_R12_INTERMED 0

#define COMPUTE_AB_BLOCK_ONLY 0
#define COMPUTE_MA_BLOCK_ONLY 0

void
R12IntEval::obs_contrib_to_VXB_gebc_vbseqobs_()
{
  if (evaluated_)
    return;
  LinearR12::ABSMethod abs_method = r12info_->abs_method();
  Ref<MessageGrp> msg = r12info_->msg();
  Ref<MemoryGrp> mem = r12info_->mem();
  Ref<ThreadGrp> thr = r12info_->thr();
  const int num_te_types = 3;
  enum te_types {eri=0, r12=1, r12t1=2};

  tim_enter("mp2-r12a intermeds");

  int me = msg->me();
  int nproc = msg->n();
  
  ExEnv::out0() << endl << indent
                << "Entered OBS A (GEBC) intermediates evaluator" << endl;
  ExEnv::out0() << incindent;

  // Do the AO->MO transform
  Ref<TwoBodyMOIntsTransform> ipjq_tform = get_tform_("(ip|jq)");
  Ref<R12IntsAcc> ipjq_acc = ipjq_tform->ints_acc();
  if (!ipjq_acc->is_committed()) {
    ipjq_tform->set_num_te_types(num_te_types);
    ipjq_tform->compute();
  }
  if (num_te_types != ipjq_acc->num_te_types())
    throw std::runtime_error("R12IntEval::obs_contrib_to_VXB_gebc() -- number of MO integral types is wrong");

  int nocc = r12info_->nocc();
  int nocc_act = r12info_->nocc_act();
  int nfzc = r12info_->nfzc();
  int nfzv = r12info_->nfzv();
  int noso = r12info_->mo_space()->rank();
  int nvir  = noso - nocc;

  /*--------------------------------
    Compute MP2-R12/A intermediates
    and collect on node0
   --------------------------------*/
  ExEnv::out0() << indent << "Begin computation of intermediates" << endl;
  tim_enter("intermediates");
  SpatialMOPairIter_eq ij_iter(r12info_->act_occ_space());
  SpatialMOPairIter_eq kl_iter(r12info_->act_occ_space());
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
  int nproc_with_ints = tasks_with_ints_(ipjq_acc,proc_with_ints);

  
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

  bool two_basis_form = (r12info_->basis() != r12info_->basis_ri());
  double pfac_xy_1, pfac_xy_2;
  if (two_basis_form &&
      ( abs_method == LinearR12::ABS_ABS ||
        abs_method == LinearR12::ABS_ABSPlus ) ) {
    pfac_xy_1 = 0.5;
    pfac_xy_2 = -0.5;
  }
  else {
    pfac_xy_1 = 0.5;
    pfac_xy_2 = 0.5;
  }
  
  if (ipjq_acc->has_access(me)) {

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
      double *klxy_buf_eri = ipjq_acc->retrieve_pair_block(k,l,R12IntsAcc::eri);
      double *klxy_buf_r12 = ipjq_acc->retrieve_pair_block(k,l,R12IntsAcc::r12);
      double *klxy_buf_r12t1 = ipjq_acc->retrieve_pair_block(k,l,R12IntsAcc::r12t1);
      double *lkxy_buf_r12t1 = ipjq_acc->retrieve_pair_block(l,k,R12IntsAcc::r12t1);
      tim_exit("MO ints retrieve");

      if (debug_)
        ExEnv::outn() << indent << "task " << me << ": obtained kl blocks" << endl;

      // Compute MP2 energies
      RefDiagSCMatrix act_occ_evals = r12info_->act_occ_space()->evals();
      RefDiagSCMatrix all_evals = r12info_->obs_space()->evals();
      double emp2_aa = 0.0;
      double emp2_ab = 0.0;
      for(int a=nocc; a<noso; a++) {
        for(int b=nocc; b<noso; b++) {
          const int ab_offset = a*noso+b;
          const int ba_offset = b*noso+a;
          const double oo_delta_ijab = 1.0/(act_occ_evals(k)+act_occ_evals(l)-all_evals(a)-all_evals(b));
          const double eri_kalb = klxy_buf_eri[ab_offset];
          const double eri_kbla = klxy_buf_eri[ba_offset];
          emp2_ab += 0.5*(eri_kalb * eri_kalb + eri_kbla * eri_kbla) * oo_delta_ijab;
          if (kl_aa != -1) {
            emp2_aa += (eri_kalb - eri_kbla) * (eri_kalb - eri_kbla) * oo_delta_ijab;
          }
        }
      }
      emp2pair_ab_.accumulate_element(kl_ab,emp2_ab);
      if (kl_ab != lk_ab)
        emp2pair_ab_.accumulate_element(lk_ab,emp2_ab);
      if (kl_aa != -1)
        emp2pair_aa_.accumulate_element(kl_aa,emp2_aa);

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
        double *ijxy_buf_r12 = ipjq_acc->retrieve_pair_block(i,j,R12IntsAcc::r12);
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

        for(int y=0;y<noso;y++) {
          double pfac_xy;
          if (y >= nocc)
            pfac_xy = pfac_xy_1;
          else
            pfac_xy = pfac_xy_2;
          for(int x=0;x<noso;x++) {

#if COMPUTE_AB_BLOCK_ONLY
            if (y < nocc || x < nocc) {
              pfac_xy = 0.0;
            }
            else {
              if (y >= nocc)
                pfac_xy = pfac_xy_1;
              else
                pfac_xy = pfac_xy_2;
            }
#endif
#if COMPUTE_MA_BLOCK_ONLY
            if ((y < nocc && x < nocc) ||
                (y >= nocc && x >= nocc)) {
              pfac_xy = 0.0;
            }
            else {
              pfac_xy = 0.5;
            }
#endif
            
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
        ipjq_acc->release_pair_block(i,j,R12IntsAcc::r12);
      }
      ipjq_acc->release_pair_block(k,l,R12IntsAcc::eri);
      ipjq_acc->release_pair_block(k,l,R12IntsAcc::r12);
      ipjq_acc->release_pair_block(k,l,R12IntsAcc::r12t1);
      ipjq_acc->release_pair_block(l,k,R12IntsAcc::r12t1);
    }
  }
  // Tasks that don't do any work here still need to create these timers
  tim_enter("MO ints retrieve");
  tim_exit("MO ints retrieve");
  tim_enter("MO ints contraction");
  tim_exit("MO ints contraction");

  tim_exit("intermediates");
  ExEnv::out0() << indent << "End of computation of intermediates" << endl;
  ipjq_acc->deactivate();
  
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
  ExEnv::out0() << indent << "Exited OBS A (GEBC) intermediates evaluator" << endl;

  tim_exit("mp2-r12a intermeds");
  checkpoint_();
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
