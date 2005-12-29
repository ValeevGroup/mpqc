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

void
R12IntEval::compute_B_gbc_1_()
{
  if (abs_method_ == LinearR12::ABS_ABS || abs_method_ == LinearR12::ABS_ABSPlus)
    throw std::runtime_error("R12IntEval::compute_B_gbc_1_() -- B(GBC1) term can only be computed using a CABS (or CABS+) approach");

  if (evaluated_)
    return;

  Ref<TwoBodyMOIntsTransform> ipjq_tform = get_tform_("(ip|jq)");
  Ref<R12IntsAcc> ijpq_acc = ipjq_tform->ints_acc();
  if (!ijpq_acc->is_committed())
    ipjq_tform->compute(intparams_);
  if (!ijpq_acc->is_active())
    ijpq_acc->activate();

  tim_enter("B(GBC1) intermediate");

  Ref<MessageGrp> msg = r12info()->msg();
  int me = msg->me();
  int nproc = msg->n();
  ExEnv::out0() << endl << indent
    << "Entered B(GBC1) intermediate evaluator" << endl;
  ExEnv::out0() << incindent;

  RefSCMatrix B_gbc1_aa = Baa_.clone();  B_gbc1_aa.assign(0.0);
  RefSCMatrix B_gbc1_ab = Bab_.clone();  B_gbc1_ab.assign(0.0);

  const Ref<MOIndexSpace>& obs_space = r12info_->refinfo()->orbs();
  const Ref<MOIndexSpace>& occ_space = r12info_->refinfo()->docc();
  const Ref<MOIndexSpace>& act_occ_space = r12info_->refinfo()->docc_act();
  Ref<MOIndexSpace> vir_space = r12info_->vir();
  Ref<MOIndexSpace> ribs_space = r12info_->ribs_space();
  form_focc_space_();
  Ref<MOIndexSpace> focc_space = focc_space_;

  const int noso = obs_space->rank();
  const int nocc = occ_space->rank();
  const int nvir = noso - nocc;
  const int nribs = ribs_space->rank();

  //
  // Do the AO->MO transform for (act_occ occ|r12|act_occ ribs) and (act_occ focc|r12|act_occ ribs)
  //
  Ref<MOIntsTransformFactory> tfactory = r12info_->tfactory();

  tfactory->set_spaces(act_occ_space,occ_space,
                       act_occ_space,ribs_space);
  Ref<TwoBodyMOIntsTransform> imjA_tform = tfactory->twobody_transform_13("(im|jA)",corrfactor_->callback());
  imjA_tform->compute(intparams_);
  Ref<R12IntsAcc> ijmA_acc = imjA_tform->ints_acc();

  tfactory->set_spaces(act_occ_space,focc_space,
                       act_occ_space,ribs_space);
  Ref<TwoBodyMOIntsTransform> iMfjA_tform = tfactory->twobody_transform_13("(iMf|jA)",corrfactor_->callback());
  iMfjA_tform->compute(intparams_);
  Ref<R12IntsAcc> ijMfA_acc = iMfjA_tform->ints_acc();
  
  SpatialMOPairIter_eq ij_iter(act_occ_space);
  SpatialMOPairIter_eq kl_iter(act_occ_space);
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

#if 1
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

    // Get (|r12|) integrals
    tim_enter("MO ints retrieve");
    double *klMfA_buf_r12 = ijMfA_acc->retrieve_pair_block(k,l,corrfactor_->tbint_type_f12());
    double *lkMfA_buf_r12 = ijMfA_acc->retrieve_pair_block(l,k,corrfactor_->tbint_type_f12());
    tim_exit("MO ints retrieve");

    if (debug_)
      ExEnv::outn() << indent << "task " << me << ": obtained kl blocks" << endl;

    for(ij_iter.start(kl+1);int(ij_iter);ij_iter.next()) {

      const int ij = ij_iter.ij();
      const int i = ij_iter.i();
      const int j = ij_iter.j();
      const int ij_aa = ij_iter.ij_aa();
      const int ij_ab = ij_iter.ij_ab();
      const int ji_ab = ij_iter.ij_ba();

      if (debug_)
        ExEnv::outn() << indent << "task " << me << ": working on (i,j) = " << i << "," << j << " " << endl;

      // Get (|r12|) integrals
      tim_enter("MO ints retrieve");
      double *ijmA_buf_r12 = ijmA_acc->retrieve_pair_block(i,j,corrfactor_->tbint_type_f12());
      double *jimA_buf_r12 = ijmA_acc->retrieve_pair_block(j,i,corrfactor_->tbint_type_f12());
      tim_exit("MO ints retrieve");

      if (debug_)
        ExEnv::outn() << indent << "task " << me << ": obtained ij blocks" << endl;

      double rr_klij = 0.0;
      double rr_lkji = 0.0;
      double rr_klji = 0.0;
      double rr_lkij = 0.0;

      const int unit_stride = 1;
      rr_klij = F77_DDOT(&nbraket,klMfA_buf_r12,&unit_stride,ijmA_buf_r12,&unit_stride);
      if (kl_ab != lk_ab && ij_ab != ji_ab) {
        rr_lkji = F77_DDOT(&nbraket,lkMfA_buf_r12,&unit_stride,jimA_buf_r12,&unit_stride);
      }
      else
        rr_lkji = rr_klij;
      B_gbc1_ab.set_element(kl_ab,ij_ab,-(rr_klij+rr_lkji));
      B_gbc1_ab.set_element(lk_ab,ji_ab,-(rr_klij+rr_lkji));
      
      if (kl_ab != lk_ab)
        rr_lkij = F77_DDOT(&nbraket,lkMfA_buf_r12,&unit_stride,ijmA_buf_r12,&unit_stride);
      else
        rr_lkij = rr_klij;
      if (ij_ab != ji_ab)
        rr_klji += F77_DDOT(&nbraket,klMfA_buf_r12,&unit_stride,jimA_buf_r12,&unit_stride);
      else
        rr_klji = rr_klij;
      B_gbc1_ab.set_element(kl_ab,ji_ab,-(rr_klji+rr_lkij));
      B_gbc1_ab.set_element(lk_ab,ij_ab,-(rr_klji+rr_lkij));

      if (ij_aa != -1 && kl_aa != -1)
        B_gbc1_aa.set_element(kl_aa,ij_aa,-(rr_klij+rr_lkji-rr_klji-rr_lkij));
        
      ijmA_acc->release_pair_block(i,j,corrfactor_->tbint_type_f12());
      ijmA_acc->release_pair_block(j,i,corrfactor_->tbint_type_f12());
    }

    ijMfA_acc->release_pair_block(k,l,corrfactor_->tbint_type_f12());
    ijMfA_acc->release_pair_block(l,k,corrfactor_->tbint_type_f12());
  }
#endif

#if 1
  //
  // Do the AO->MO transform for (act_occ focc|r12|act_occ vir)
  //
  tfactory->set_spaces(act_occ_space,focc_space,
                       act_occ_space,vir_space);
  Ref<TwoBodyMOIntsTransform> iMfja_tform = tfactory->twobody_transform_13("(iMf|ja)",corrfactor_->callback());
  iMfja_tform->compute(intparams_);
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
    const int lk_ab = kl_iter.ij_ba();

    if (debug_)
      ExEnv::outn() << indent << "task " << me << ": working on (k,l) = " << k << "," << l << " " << endl;

    // Get (|r12|) integrals
    tim_enter("MO ints retrieve");
    double *klMfa_buf_r12 = ijMfa_acc->retrieve_pair_block(k,l,corrfactor_->tbint_type_f12());
    double *lkMfa_buf_r12 = ijMfa_acc->retrieve_pair_block(l,k,corrfactor_->tbint_type_f12());
    tim_exit("MO ints retrieve");

    if (debug_)
      ExEnv::outn() << indent << "task " << me << ": obtained kl blocks" << endl;

    for(ij_iter.start(kl+1);int(ij_iter);ij_iter.next()) {

      const int ij = ij_iter.ij();
      const int i = ij_iter.i();
      const int j = ij_iter.j();
      const int ij_aa = ij_iter.ij_aa();
      const int ij_ab = ij_iter.ij_ab();
      const int ji_ab = ij_iter.ij_ba();

      if (debug_)
        ExEnv::outn() << indent << "task " << me << ": working on (i,j) = " << i << "," << j << " " << endl;

      // Get (|r12|) integrals
      tim_enter("MO ints retrieve");
      double *ijpq_buf_r12 = ijpq_acc->retrieve_pair_block(i,j,corrfactor_->tbint_type_f12());
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
          rr_lkji = lkMfa_buf_r12[ma]*ijpq_buf_r12[am_offset];
          B_gbc1_ab.accumulate_element(kl_ab,ij_ab,-(rr_klij+rr_lkji));
          if (kl_ab != lk_ab && ij_ab != ji_ab) {
            B_gbc1_ab.accumulate_element(lk_ab,ji_ab,-(rr_lkji+rr_klij));
          }
          
          rr_lkij = lkMfa_buf_r12[ma]*ijpq_buf_r12[ma_offset];
          rr_klji = klMfa_buf_r12[ma]*ijpq_buf_r12[am_offset];
          if (kl_ab != lk_ab)
            B_gbc1_ab.accumulate_element(lk_ab,ij_ab,-(rr_lkij+rr_klji));
          if (ij_ab != ji_ab)
            B_gbc1_ab.accumulate_element(kl_ab,ji_ab,-(rr_klji+rr_lkij));
          
          if (ij_aa != -1 && kl_aa != -1)
            B_gbc1_aa.accumulate_element(kl_aa,ij_aa,-(rr_klij-rr_lkij-rr_klji+rr_lkji));
        }
      
      ijpq_acc->release_pair_block(i,j,corrfactor_->tbint_type_f12());
    }

    ijMfa_acc->release_pair_block(k,l,corrfactor_->tbint_type_f12());
    ijMfa_acc->release_pair_block(l,k,corrfactor_->tbint_type_f12());
  }
#endif

  if (debug_ > 1) {
    B_gbc1_aa.print("Alpha-alpha B(GBC1) contribution");
    B_gbc1_ab.print("Alpha-beta B(GBC1) contribution");
  }
  // Symmetrize the B contribution
  B_gbc1_aa.scale(0.5);
  B_gbc1_ab.scale(0.5);
  RefSCMatrix B_gbc1_aa_t = B_gbc1_aa.t();
  Baa_.accumulate(B_gbc1_aa); Baa_.accumulate(B_gbc1_aa_t);
  RefSCMatrix B_gbc1_ab_t = B_gbc1_ab.t();
  Bab_.accumulate(B_gbc1_ab); Bab_.accumulate(B_gbc1_ab_t);

  globally_sum_intermeds_();

  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Exited B(GBC1) intermediate evaluator" << endl;

  tim_exit("B(GBC1) intermediate");
}


void
R12IntEval::compute_B_gbc_2_()
{
  if (abs_method_ == LinearR12::ABS_ABS || abs_method_ == LinearR12::ABS_ABSPlus)
    throw std::runtime_error("R12IntEval::compute_B_gbc_2_() -- B(GBC2) term can only be computed using a CABS (or CABS+) approach");

  if (evaluated_)
    return;

  Ref<TwoBodyMOIntsTransform> ipjq_tform = get_tform_("(ip|jq)");
  Ref<R12IntsAcc> ijpq_acc = ipjq_tform->ints_acc();
  if (!ijpq_acc->is_committed())
    ipjq_tform->compute(intparams_);
  if (!ijpq_acc->is_active())
    ijpq_acc->activate();

  tim_enter("B(GBC2) intermediate");

  Ref<MessageGrp> msg = r12info_->msg();
  int me = msg->me();
  int nproc = msg->n();
  ExEnv::out0() << endl << indent
    << "Entered B(GBC2) intermediate evaluator" << endl;
  ExEnv::out0() << incindent;

  RefSCMatrix X_ijklF_ab = Bab_.clone();
  RefSCMatrix B_gbc2_aa = Baa_.clone();
  RefSCMatrix B_gbc2_ab = Bab_.clone();
  X_ijklF_ab.assign(0.0);
  B_gbc2_aa.assign(0.0);
  B_gbc2_ab.assign(0.0);
  
  const Ref<MOIndexSpace>& obs_space = r12info_->refinfo()->orbs();
  const Ref<MOIndexSpace>& occ_space = r12info_->refinfo()->docc();
  const Ref<MOIndexSpace>& act_occ_space = r12info_->refinfo()->docc_act();
  Ref<MOIndexSpace> ribs_space = r12info_->ribs_space();
  form_factocc_space_();
  Ref<MOIndexSpace> factocc_space = factocc_space_;

  const int nocc = occ_space->rank();
  const int nribs = ribs_space->rank();

  // compute r_{12}^2 operator in act.occ.pair/act.occ.-focc. basis
  RefSCMatrix R2 = compute_r2_(act_occ_space,act_occ_space,act_occ_space,factocc_space);
  // Compute contribution X += (r^2)_{ij}^{k l_f}
  if (me == 0)
    X_ijklF_ab.accumulate(R2);

  //
  // Compute contribution X -= r_{ij}^{\alpha'm} r_{m\alpha'}^{k l_f}
  //                         + r_{ji}^{\alpha'm} r_{\alpha'm}^{k l_f}
  //
  Ref<MOIntsTransformFactory> tfactory = r12info_->tfactory();

  tfactory->set_spaces(act_occ_space,occ_space,
                       act_occ_space,ribs_space);
  Ref<TwoBodyMOIntsTransform> imjA_tform = tfactory->twobody_transform_13("(im|jA)",corrfactor_->callback());
  imjA_tform->compute(intparams_);
  Ref<R12IntsAcc> ijmA_acc = imjA_tform->ints_acc();

  tfactory->set_spaces(act_occ_space,occ_space,
                       factocc_space,ribs_space);
  Ref<TwoBodyMOIntsTransform> kmlfA_tform = tfactory->twobody_transform_13("(km|lfA)",corrfactor_->callback());
  kmlfA_tform->compute(intparams_);
  Ref<R12IntsAcc> klfmA_acc = kmlfA_tform->ints_acc();

  tfactory->set_spaces(factocc_space,occ_space,
                       act_occ_space,ribs_space);
  Ref<TwoBodyMOIntsTransform> lfmkA_tform = tfactory->twobody_transform_13("(lfm|kA)",corrfactor_->callback());
  lfmkA_tform->compute(intparams_);
  Ref<R12IntsAcc> lfkmA_acc = lfmkA_tform->ints_acc();

  // Compute the number of tasks that have full access to the integrals
  // and split the work among them
  vector<int> proc_with_ints;
  int nproc_with_ints = tasks_with_ints_(ijmA_acc,proc_with_ints);

  SpatialMOPairIter_eq ij_iter(act_occ_space);
  SpatialMOPairIter_eq kl_iter(act_occ_space);

  int nbraket = nocc*nribs;
#if 1
  for(kl_iter.start();int(kl_iter);kl_iter.next()) {

    const int kl = kl_iter.ij();
    // Figure out if this task will handle this kl
    int kl_proc = kl%nproc_with_ints;
    if (kl_proc != proc_with_ints[me])
      continue;
    const int k = kl_iter.i();
    const int l = kl_iter.j();
    const int kl_ab = kl_iter.ij_ab();
    const int lk_ab = kl_iter.ij_ba();

    if (debug_)
      ExEnv::outn() << indent << "task " << me << ": working on (k,l) = " << k << "," << l << " " << endl;

    // Get (|r12|) integrals
    tim_enter("MO ints retrieve");
    double *klfmA_buf_r12 = klfmA_acc->retrieve_pair_block(k,l,corrfactor_->tbint_type_f12());
    double *lfkmA_buf_r12 = lfkmA_acc->retrieve_pair_block(l,k,corrfactor_->tbint_type_f12());
    double *lkfmA_buf_r12 = klfmA_acc->retrieve_pair_block(l,k,corrfactor_->tbint_type_f12());
    double *kflmA_buf_r12 = lfkmA_acc->retrieve_pair_block(k,l,corrfactor_->tbint_type_f12());
    tim_exit("MO ints retrieve");

    if (debug_)
      ExEnv::outn() << indent << "task " << me << ": obtained kl blocks" << endl;

    for(ij_iter.start();int(ij_iter);ij_iter.next()) {

      const int ij = ij_iter.ij();
      const int i = ij_iter.i();
      const int j = ij_iter.j();
      const int ij_ab = ij_iter.ij_ab();
      const int ji_ab = ij_iter.ij_ba();

      if (debug_)
        ExEnv::outn() << indent << "task " << me << ": working on (i,j) = " << i << "," << j << " " << endl;

      // Get (|r12|) integrals
      tim_enter("MO ints retrieve");
      double *ijmA_buf_r12 = ijmA_acc->retrieve_pair_block(i,j,corrfactor_->tbint_type_f12());
      double *jimA_buf_r12 = ijmA_acc->retrieve_pair_block(j,i,corrfactor_->tbint_type_f12());
      tim_exit("MO ints retrieve");

      if (debug_)
        ExEnv::outn() << indent << "task " << me << ": obtained ij blocks" << endl;

      double X_ijklf = 0.0;
      double X_jiklf = 0.0;
      double X_ijlkf = 0.0;
      double X_jilkf = 0.0;

      const int unit_stride = 1;
      X_ijklf += F77_DDOT(&nbraket,lfkmA_buf_r12,&unit_stride,jimA_buf_r12,&unit_stride);
      X_ijklf += F77_DDOT(&nbraket,klfmA_buf_r12,&unit_stride,ijmA_buf_r12,&unit_stride);
      X_ijklF_ab.accumulate_element(ij_ab,kl_ab,-X_ijklf);
      if (kl_ab != lk_ab) {
        X_ijlkf += F77_DDOT(&nbraket,kflmA_buf_r12,&unit_stride,jimA_buf_r12,&unit_stride);
        X_ijlkf += F77_DDOT(&nbraket,lkfmA_buf_r12,&unit_stride,ijmA_buf_r12,&unit_stride);
        X_ijklF_ab.accumulate_element(ij_ab,lk_ab,-X_ijlkf);
      }
      if (ij_ab != ji_ab) {
        X_jiklf += F77_DDOT(&nbraket,lfkmA_buf_r12,&unit_stride,ijmA_buf_r12,&unit_stride);
        X_jiklf += F77_DDOT(&nbraket,klfmA_buf_r12,&unit_stride,jimA_buf_r12,&unit_stride);
        X_ijklF_ab.accumulate_element(ji_ab,kl_ab,-X_jiklf);
        if (kl_ab != lk_ab) {
          X_jilkf += F77_DDOT(&nbraket,kflmA_buf_r12,&unit_stride,ijmA_buf_r12,&unit_stride);
          X_jilkf += F77_DDOT(&nbraket,lkfmA_buf_r12,&unit_stride,jimA_buf_r12,&unit_stride);
          X_ijklF_ab.accumulate_element(ji_ab,lk_ab,-X_jilkf);
        }
      }
      
      ijmA_acc->release_pair_block(i,j,corrfactor_->tbint_type_f12());
      ijmA_acc->release_pair_block(j,i,corrfactor_->tbint_type_f12());
    }

    klfmA_acc->release_pair_block(k,l,corrfactor_->tbint_type_f12());
    lfkmA_acc->release_pair_block(l,k,corrfactor_->tbint_type_f12());
    klfmA_acc->release_pair_block(l,k,corrfactor_->tbint_type_f12());
    lfkmA_acc->release_pair_block(k,l,corrfactor_->tbint_type_f12());
  }
#endif

#if 1
  //
  // Compute contribution X -= r_{ij}^{pq} r_{pq}^{k l_f}
  //
  tfactory->set_spaces(act_occ_space,obs_space,
                       factocc_space,obs_space);
  Ref<TwoBodyMOIntsTransform> kplfq_tform = tfactory->twobody_transform_13("(kp|lfq)",corrfactor_->callback());
  kplfq_tform->compute(intparams_);
  Ref<R12IntsAcc> klfpq_acc = kplfq_tform->ints_acc();

  nbraket = obs_space->rank() * obs_space->rank();
  for(kl_iter.start();int(kl_iter);kl_iter.next()) {

    const int kl = kl_iter.ij();
    // Figure out if this task will handle this kl
    int kl_proc = kl%nproc_with_ints;
    if (kl_proc != proc_with_ints[me])
      continue;
    const int k = kl_iter.i();
    const int l = kl_iter.j();
    const int kl_ab = kl_iter.ij_ab();
    const int lk_ab = kl_iter.ij_ba();

    if (debug_)
      ExEnv::outn() << indent << "task " << me << ": working on (k,l) = " << k << "," << l << " " << endl;

    // Get (|r12|) integrals
    tim_enter("MO ints retrieve");
    double *klfpq_buf_r12 = klfpq_acc->retrieve_pair_block(k,l,corrfactor_->tbint_type_f12());
    double *lkfpq_buf_r12 = klfpq_acc->retrieve_pair_block(l,k,corrfactor_->tbint_type_f12());
    tim_exit("MO ints retrieve");

    if (debug_)
      ExEnv::outn() << indent << "task " << me << ": obtained kl blocks" << endl;

    for(ij_iter.start();int(ij_iter);ij_iter.next()) {

      const int ij = ij_iter.ij();
      const int i = ij_iter.i();
      const int j = ij_iter.j();
      const int ij_ab = ij_iter.ij_ab();
      const int ji_ab = ij_iter.ij_ba();

      if (debug_)
        ExEnv::outn() << indent << "task " << me << ": working on (i,j) = " << i << "," << j << " " << endl;

      // Get (|r12|) integrals
      tim_enter("MO ints retrieve");
      double *ijpq_buf_r12 = ijpq_acc->retrieve_pair_block(i,j,corrfactor_->tbint_type_f12());
      double *jipq_buf_r12 = ijpq_acc->retrieve_pair_block(j,i,corrfactor_->tbint_type_f12());
      tim_exit("MO ints retrieve");

      if (debug_)
        ExEnv::outn() << indent << "task " << me << ": obtained ij blocks" << endl;

      double X_ijklf = 0.0;
      double X_jiklf = 0.0;
      double X_ijlkf = 0.0;
      double X_jilkf = 0.0;

      const int unit_stride = 1;
      X_ijklf += F77_DDOT(&nbraket,klfpq_buf_r12,&unit_stride,ijpq_buf_r12,&unit_stride);
      X_ijklF_ab.accumulate_element(ij_ab,kl_ab,-X_ijklf);
      if (kl_ab != lk_ab) {
        X_ijlkf += F77_DDOT(&nbraket,lkfpq_buf_r12,&unit_stride,ijpq_buf_r12,&unit_stride);
        X_ijklF_ab.accumulate_element(ij_ab,lk_ab,-X_ijlkf);
      }
      if (ij_ab != ji_ab) {
        X_jiklf += F77_DDOT(&nbraket,klfpq_buf_r12,&unit_stride,jipq_buf_r12,&unit_stride);
        X_ijklF_ab.accumulate_element(ji_ab,kl_ab,-X_jiklf);
        if (kl_ab != lk_ab) {
          X_jilkf += F77_DDOT(&nbraket,lkfpq_buf_r12,&unit_stride,jipq_buf_r12,&unit_stride);
          X_ijklF_ab.accumulate_element(ji_ab,lk_ab,-X_jilkf);
        }
      }
      
      ijpq_acc->release_pair_block(i,j,corrfactor_->tbint_type_f12());
      ijpq_acc->release_pair_block(j,i,corrfactor_->tbint_type_f12());
    }

    klfpq_acc->release_pair_block(k,l,corrfactor_->tbint_type_f12());
    klfpq_acc->release_pair_block(l,k,corrfactor_->tbint_type_f12());
  }
  globally_sum_scmatrix_(X_ijklF_ab);
#endif

  //
  // Compute B_gbc2 = X_ijklF + X_jilkF :
  // B_gbc2_ab_ijkl = X_ijklF_ab + X_jilkF_ab
  // B_gbc2_aa_ijkl = X_ijklF_aa + X_jilkF_aa = X_ijklF_ab - X_jiklF_ab + X_jilkF_ab - X_ijlkF_ab
  //

  for(kl_iter.start();int(kl_iter);kl_iter.next()) {

    const int kl_aa = kl_iter.ij_aa();
    const int kl_ab = kl_iter.ij_ab();
    const int lk_ab = kl_iter.ij_ba();

    for(ij_iter.start();int(ij_iter);ij_iter.next()) {

      const int ij_aa = ij_iter.ij_aa();
      const int ij_ab = ij_iter.ij_ab();
      const int ji_ab = ij_iter.ij_ba();

      const double B_ab_ijkl = X_ijklF_ab.get_element(ij_ab,kl_ab) + X_ijklF_ab.get_element(ji_ab,lk_ab);
      const double B_ab_ijlk = X_ijklF_ab.get_element(ij_ab,lk_ab) + X_ijklF_ab.get_element(ji_ab,kl_ab);
      const double B_ab_jikl = B_ab_ijlk;
      const double B_ab_jilk = B_ab_ijkl;

      B_gbc2_ab.set_element( ij_ab, kl_ab, B_ab_ijkl);
      if (kl_ab != lk_ab)
        B_gbc2_ab.set_element( ij_ab, lk_ab, B_ab_ijlk);
      if (ij_ab != ji_ab) {
        B_gbc2_ab.set_element( ji_ab, kl_ab, B_ab_jikl);
        if (kl_ab != lk_ab)
          B_gbc2_ab.set_element( ji_ab, lk_ab, B_ab_jilk);
      }

      if (ij_aa != -1 && kl_aa != -1) {
        B_gbc2_aa.set_element( ij_aa, kl_aa, B_ab_ijkl - B_ab_jikl);
      }

    }
  }

  if (debug_ > 1) {
    B_gbc2_aa.print("Alpha-alpha B(GBC2) contribution");
    B_gbc2_ab.print("Alpha-beta B(GBC2) contribution");
  }
  // Symmetrize the B contribution
  B_gbc2_aa.scale(0.5);
  B_gbc2_ab.scale(0.5);
  RefSCMatrix B_gbc2_aa_t = B_gbc2_aa.t();
  Baa_.accumulate(B_gbc2_aa); Baa_.accumulate(B_gbc2_aa_t);
  RefSCMatrix B_gbc2_ab_t = B_gbc2_ab.t();
  Bab_.accumulate(B_gbc2_ab); Bab_.accumulate(B_gbc2_ab_t);

  globally_sum_intermeds_();

  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Exited B(GBC2) intermediate evaluator" << endl;

  tim_exit("B(GBC2) intermediate");
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
