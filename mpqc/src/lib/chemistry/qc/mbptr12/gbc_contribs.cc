//
// gbc_contribs.cc
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
#include <string.h>
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
#include <chemistry/qc/mbptr12/blas.h>
#include <chemistry/qc/mbptr12/r12ia.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>
#include <chemistry/qc/mbptr12/pairiter.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>

using namespace std;
using namespace sc;

// Computes the "extra" contribution to B that appears during the initial rearrangement

void
R12IntEval::compute_B_gbc_1_()
{
  if (evaluated_)
    return;

  Ref<R12IntsAcc> ijpq_acc = ipjq_tform_->ints_acc();

  if (!ijpq_acc->is_committed())
    throw std::runtime_error("R12IntEval::compute_B_gbc_1_() -- ipjq_tform_ hasn't been computed yet");
  if (!ijpq_acc->is_active())
    ijpq_acc->activate();
  
  tim_enter("B(GBC1) intermediate");

  const int num_te_types = 2;
  
  Ref<MessageGrp> msg = r12info()->msg();
  int me = msg->me();
  int nproc = msg->n();
  ExEnv::out0() << endl << indent
    << "Entered B(GBC1) intermediate evaluator" << endl;
  ExEnv::out0() << indent << scprintf("nproc = %i", nproc) << endl;

  RefSCMatrix B_gbc1_aa = Baa_.clone();  B_gbc1_aa.assign(0.0);
  RefSCMatrix B_gbc1_ab = Bab_.clone();  B_gbc1_ab.assign(0.0);

  Ref<MOIndexSpace> mo_space = r12info_->obs_space();
  Ref<MOIndexSpace> occ_space = r12info_->occ_space();
  Ref<MOIndexSpace> act_occ_space = r12info_->act_occ_space();
  Ref<MOIndexSpace> vir_space = r12info_->vir_space();

  // compute the Fock matrix between the complement and all occupieds and
  // create the new Fock-weighted space
  Ref<MOIndexSpace> ribs_space = r12info_->ribs_space();
  RefSCMatrix F_ri_o = fock_(r12info_->occ_space(),ribs_space,occ_space);
  F_ri_o.print("Fock matrix (RI-BS/occ.)");
  Ref<MOIndexSpace> focc_space = new MOIndexSpace("Fock-weighted occupied MOs sorted by energy",
                                                  occ_space, ribs_space->coefs()*F_ri_o, ribs_space->basis());

  const int noso = r12info()->noso();
  const int nocc = r12info_->nocc();
  const int nvir = noso - nocc;
  const int nribs = ribs_space->rank();

  //
  // Do the AO->MO transform for (act_occ occ|r12|act_occ ribs) and (act_occ focc|r12|act_occ ribs)
  //
  Ref<MOIntsTransformFactory> tfactory = r12info_->tfactory();

  tfactory->set_spaces(act_occ_space,occ_space,
                       act_occ_space,ribs_space);
  Ref<TwoBodyMOIntsTransform> imjA_tform = tfactory->twobody_transform_13("(im|jA)");
  imjA_tform->set_num_te_types(num_te_types);
  imjA_tform->compute();
  Ref<R12IntsAcc> ijmA_acc = imjA_tform->ints_acc();

  tfactory->set_spaces(act_occ_space,focc_space,
                       act_occ_space,ribs_space);
  Ref<TwoBodyMOIntsTransform> iMfjA_tform = tfactory->twobody_transform_13("(iMf|jA)");
  iMfjA_tform->set_num_te_types(num_te_types);
  iMfjA_tform->compute();
  Ref<R12IntsAcc> ijMfA_acc = iMfjA_tform->ints_acc();
  
  MOPairIter_SD ij_iter(act_occ_space);
  MOPairIter_SD kl_iter(act_occ_space);
  int naa = ij_iter.nij_aa();          // Number of alpha-alpha pairs (i > j)
  int nab = ij_iter.nij_ab();          // Number of alpha-beta pairs
  if (debug_) {
    ExEnv::out0() << indent << "naa = " << naa << endl;
    ExEnv::out0() << indent << "nab = " << nab << endl;
  }

  // Compute the number of tasks that have full access to the integrals
  // and split the work among them
  vector<int> proc_with_ints;
  int nproc_with_ints = tasks_with_ints_(ijMfA_acc,proc_with_ints);

  // Compute the first half of the term
  const int nbraket = nocc*nribs;

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
    const int lk_ab = kl_iter.ji_ab();

    if (debug_)
      ExEnv::outn() << indent << "task " << me << ": working on (k,l) = " << k << "," << l << " " << endl;

    // Get (|r12|) integrals
    tim_enter("MO ints retrieve");
    double *klMfA_buf_r12 = ijMfA_acc->retrieve_pair_block(k,l,R12IntsAcc::r12);
    double *lkMfA_buf_r12 = ijMfA_acc->retrieve_pair_block(l,k,R12IntsAcc::r12);
    tim_exit("MO ints retrieve");

    if (debug_)
      ExEnv::outn() << indent << "task " << me << ": obtained kl blocks" << endl;

    for(ij_iter.start(kl+1);int(ij_iter);ij_iter.next()) {

      const int ij = ij_iter.ij();
      const int i = ij_iter.i();
      const int j = ij_iter.j();
      const int ij_aa = ij_iter.ij_aa();
      const int ij_ab = ij_iter.ij_ab();
      const int ji_ab = ij_iter.ji_ab();

      if (debug_)
        ExEnv::outn() << indent << "task " << me << ": working on (i,j) = " << i << "," << j << " " << endl;

      // Get (|r12|) integrals
      tim_enter("MO ints retrieve");
      double *ijmA_buf_r12 = ijmA_acc->retrieve_pair_block(i,j,R12IntsAcc::r12);
      double *jimA_buf_r12 = ijmA_acc->retrieve_pair_block(j,i,R12IntsAcc::r12);
      tim_exit("MO ints retrieve");

      if (debug_)
        ExEnv::outn() << indent << "task " << me << ": obtained ij blocks" << endl;

      double rr_klij = 0.0;
      double rr_lkij = 0.0;
      double rr_klji = 0.0;
      double rr_lkji = 0.0;

      const int unit_stride = 1;
      rr_klij = F77_DDOT(&nbraket,klMfA_buf_r12,&unit_stride,ijmA_buf_r12,&unit_stride);
      B_gbc1_ab.set_element(kl_ab,ij_ab,-rr_klij);
      if (kl_ab == lk_ab) {
        rr_lkij = F77_DDOT(&nbraket,lkMfA_buf_r12,&unit_stride,ijmA_buf_r12,&unit_stride);
        B_gbc1_ab.set_element(lk_ab,ij_ab,-rr_lkij);
      }
      if (ij_ab == ji_ab) {
        rr_klji = F77_DDOT(&nbraket,klMfA_buf_r12,&unit_stride,jimA_buf_r12,&unit_stride);
        B_gbc1_ab.set_element(kl_ab,ji_ab,-rr_klji);
      }
      if (kl_ab == lk_ab && ij_ab == ji_ab) {
        rr_lkji = F77_DDOT(&nbraket,lkMfA_buf_r12,&unit_stride,jimA_buf_r12,&unit_stride);
        B_gbc1_ab.set_element(lk_ab,ji_ab,-rr_lkji);
      }

      if (ij_aa != -1 && kl_aa != -1)
        B_gbc1_aa.set_element(kl_aa,ij_aa,-(rr_klij-rr_lkij-rr_klji+rr_lkji));
        
      ijmA_acc->release_pair_block(i,j,R12IntsAcc::r12);
      ijmA_acc->release_pair_block(j,i,R12IntsAcc::r12);
    }

    ijMfA_acc->release_pair_block(k,l,R12IntsAcc::r12);
    ijMfA_acc->release_pair_block(l,k,R12IntsAcc::r12);
  }

  //
  // Do the AO->MO transform for (act_occ focc|r12|act_occ vir)
  //
  tfactory->set_spaces(act_occ_space,focc_space,
                       act_occ_space,vir_space);
  Ref<TwoBodyMOIntsTransform> iMfja_tform = tfactory->twobody_transform_13("(iMf|ja)");
  iMfja_tform->set_num_te_types(num_te_types);
  iMfja_tform->compute();
  Ref<R12IntsAcc> ijMfa_acc = iMfja_tform->ints_acc();

  nproc_with_ints = tasks_with_ints_(ijMfa_acc,proc_with_ints);

  // Compute the second half of the term
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
    const int lk_ab = kl_iter.ji_ab();

    if (debug_)
      ExEnv::outn() << indent << "task " << me << ": working on (k,l) = " << k << "," << l << " " << endl;

    // Get (|r12|) integrals
    tim_enter("MO ints retrieve");
    double *klMfa_buf_r12 = ijMfa_acc->retrieve_pair_block(k,l,R12IntsAcc::r12);
    double *lkMfa_buf_r12 = ijMfa_acc->retrieve_pair_block(l,k,R12IntsAcc::r12);
    tim_exit("MO ints retrieve");

    if (debug_)
      ExEnv::outn() << indent << "task " << me << ": obtained kl blocks" << endl;

    for(ij_iter.start(kl+1);int(ij_iter);ij_iter.next()) {

      const int ij = ij_iter.ij();
      const int i = ij_iter.i();
      const int j = ij_iter.j();
      const int ij_aa = ij_iter.ij_aa();
      const int ij_ab = ij_iter.ij_ab();
      const int ji_ab = ij_iter.ji_ab();

      if (debug_)
        ExEnv::outn() << indent << "task " << me << ": working on (i,j) = " << i << "," << j << " " << endl;

      // Get (|r12|) integrals
      tim_enter("MO ints retrieve");
      double *ijpq_buf_r12 = ijpq_acc->retrieve_pair_block(i,j,R12IntsAcc::r12);
      tim_exit("MO ints retrieve");

      if (debug_)
        ExEnv::outn() << indent << "task " << me << ": obtained ij blocks" << endl;

      int ma = 0;
      for(int m=0; m<nocc; m++)
        for(int a=0; a<nvir; a++, ma++) {
          
          const int aa = a+nocc;
          const int ma_offset = m*noso+aa;
          const int am_offset = aa*noso+m;

          double rr_klij = 0.0;
          double rr_lkij = 0.0;
          double rr_klji = 0.0;
          double rr_lkji = 0.0;
          
          rr_klij = klMfa_buf_r12[ma]*ijpq_buf_r12[ma_offset];
          B_gbc1_ab.accumulate_element(kl_ab,ij_ab,-rr_klij);
          
          if (kl_ab == lk_ab) {
            rr_lkij = lkMfa_buf_r12[ma]*ijpq_buf_r12[ma_offset];
            B_gbc1_ab.accumulate_element(lk_ab,ij_ab,-rr_lkij);
          }
          if (ij_ab == ji_ab) {
            rr_klji = klMfa_buf_r12[ma]*ijpq_buf_r12[am_offset];
            B_gbc1_ab.accumulate_element(kl_ab,ji_ab,-rr_klji);
          }
          if (kl_ab == lk_ab && ij_ab == ji_ab) {
            rr_lkji = lkMfa_buf_r12[ma]*ijpq_buf_r12[am_offset];
            B_gbc1_ab.accumulate_element(lk_ab,ji_ab,-rr_lkji);
          }
          
          if (ij_aa != -1 && kl_aa != -1)
            B_gbc1_aa.accumulate_element(kl_aa,ij_aa,-(rr_klij-rr_lkij-rr_klji+rr_lkji));
        }
      
      ijpq_acc->release_pair_block(i,j,R12IntsAcc::r12);
    }

    ijMfa_acc->release_pair_block(k,l,R12IntsAcc::r12);
    ijMfa_acc->release_pair_block(l,k,R12IntsAcc::r12);
  }

  ExEnv::out0() << indent << "Exited B(GBC1) intermediate evaluator" << endl;

  // Symmetrize the B contribution
  B_gbc1_aa.scale(0.5);
  B_gbc1_ab.scale(0.5);
  B_gbc1_aa.print("Alpha-alpha B(GBC1) contribution");
  B_gbc1_ab.print("Alpha-beta B(GBC1) contribution");
  Baa_.accumulate(B_gbc1_aa); Baa_.accumulate(B_gbc1_aa.t());
  Bab_.accumulate(B_gbc1_ab); Bab_.accumulate(B_gbc1_ab.t());

  globally_sum_intermeds_();

  tim_exit("B(GBC1) intermediate");
}


////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
