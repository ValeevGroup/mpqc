//
// compute_amps.cc
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
#include <math/scmat/local.h>
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
R12IntEval::compute_T2_vbsneqobs_()
{
  Ref<TwoBodyMOIntsTransform> iajb_tform = get_tform_("(ia|jb)");
  if (iajb_tform->space1() != iajb_tform->space3())
    throw std::runtime_error("R12IntEval::compute_T2_vbsneqobs_() -- wrong type of transform is provided (space1 != space3)");
  if (iajb_tform->space2() != iajb_tform->space4())
    throw std::runtime_error("R12IntEval::compute_T2_vbsneqobs_() -- wrong type of transform is provided (space2 != space4)");
  if (iajb_tform->space1() != r12info_->refinfo()->docc_act())
    throw std::runtime_error("R12IntEval::compute_T2_vbsneqobs_() -- wrong type of transform is provided (space1 != act_occ)");
  if (iajb_tform->space2() != r12info_->vir_act())
    throw std::runtime_error("R12IntEval::compute_T2_vbsneqobs_() -- wrong type of transform is provided (space1 != act_vir)");
  Ref<R12IntsAcc> ijab_acc = iajb_tform->ints_acc();
  if (!ijab_acc->is_committed())
    iajb_tform->compute(intparams_);
  if (!ijab_acc->is_active())
    ijab_acc->activate();

  tim_enter("mp2 t2 amplitudes");

  Ref<MessageGrp> msg = r12info_->msg();
  int me = msg->me();
  int nproc = msg->n();
  ExEnv::out0() << endl << indent
	       << "Entered MP2 T2 amplitude evaluator" << endl;
  ExEnv::out0() << incindent;

  const Ref<MOIndexSpace>& act_occ_space = r12info_->refinfo()->docc_act();
  const Ref<MOIndexSpace>& act_vir_space = r12info_->vir_act();
  const int nactvir = act_vir_space->rank();
  RefDiagSCMatrix act_occ_evals = act_occ_space->evals();
  RefDiagSCMatrix act_vir_evals = act_vir_space->evals();
  
  SpatialMOPairIter_eq ij_iter(act_occ_space);
  SpatialMOPairIter_eq ab_iter(act_vir_space);
  int naa = ij_iter.nij_aa();          // Number of alpha-alpha pairs (i > j)
  int nab = ij_iter.nij_ab();          // Number of alpha-beta pairs
  if (debug_) {
    ExEnv::out0() << indent << "naa = " << naa << endl;
    ExEnv::out0() << indent << "nab = " << nab << endl;
  }
  
  // Compute the number of tasks that have full access to the integrals
  // and split the work among them
  vector<int> proc_with_ints;
  int nproc_with_ints = tasks_with_ints_(ijab_acc,proc_with_ints);
  
  //////////////////////////////////////////////////////////////
  //
  // Evaluation of the MP2 T2 amplitudes proceeds as follows:
  //
  //    loop over batches of ij,
  //      load (ijxy)=(ix|jy) into memory
  //
  //      loop over xy, 0<=x<nvir_act, 0<=y<nvir_act
  //        compute T2_aa[ij][xy] = [ (ijxy) - (ijyx) ] / denom
  //        compute T2_ab[ij][xy] = [ (ijxy) ] / denom
  //      end xy loop
  //    end ij loop
  //
  /////////////////////////////////////////////////////////////////////////////////

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

    // Get (|1/r12|) integrals
    tim_enter("MO ints retrieve");
    double *ijxy_buf_eri = ijab_acc->retrieve_pair_block(i,j,corrfactor_->tbint_type_eri());
    tim_exit("MO ints retrieve");

    if (debug_)
      ExEnv::outn() << indent << "task " << me << ": obtained ij blocks" << endl;

    // Compute MP2 energies
    double T2_aa_ijab = 0.0;
    double T2_ab_ijab = 0.0;

    for(ab_iter.start();int(ab_iter);ab_iter.next()) {

      const int a = ab_iter.i();
      const int b = ab_iter.j();
      const int ab_aa = ab_iter.ij_aa();
      const int ab_ab = ab_iter.ij_ab();
      const int ba_ab = ab_iter.ij_ba();

      const int ab_offset = a*nactvir + b;
      const int ba_offset = b*nactvir + a;
      const double oo_delta_ijab = -1.0/(-act_occ_evals(i)-act_occ_evals(j)+act_vir_evals(a)+act_vir_evals(b));
      const double eri_iajb = ijxy_buf_eri[ab_offset];
      const double eri_ibja = ijxy_buf_eri[ba_offset];
      const double T2_ab_ijab = eri_iajb * oo_delta_ijab;
      const double T2_ab_ijba = eri_ibja * oo_delta_ijab;
      T2ab_.set_element(ij_ab,ab_ab,T2_ab_ijab);
      T2ab_.set_element(ji_ab,ba_ab,T2_ab_ijab);
      T2ab_.set_element(ji_ab,ab_ab,T2_ab_ijba);
      T2ab_.set_element(ij_ab,ba_ab,T2_ab_ijba);

      if (ij_aa != -1 && ab_aa != -1) {
        const double T2_aa_ijab = (eri_iajb - eri_ibja) * oo_delta_ijab;
        T2aa_.set_element(ij_aa,ab_aa,T2_aa_ijab);
      }

    }
    
    ijab_acc->release_pair_block(i,j,corrfactor_->tbint_type_eri());
  }

  globally_sum_scmatrix_(T2aa_);
  globally_sum_scmatrix_(T2ab_);

  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Exited MP2 T2 amplitude evaluator" << endl;
  tim_exit("mp2 t2 amplitudes");
}

void
R12IntEval::compute_R_vbsneqobs_(const Ref<TwoBodyMOIntsTransform>& ipjq_tform, RefSCMatrix& Raa, RefSCMatrix& Rab)
{
  Ref<R12IntsAcc> ijpq_acc = ipjq_tform->ints_acc();
  if (ipjq_tform->space1() != ipjq_tform->space3())
    throw std::runtime_error("R12IntEval::compute_R_vbsneqobs_() -- wrong type of transform is provided (space1 != space3)");
  if (ipjq_tform->space2() != ipjq_tform->space4())
    throw std::runtime_error("R12IntEval::compute_R_vbsneqobs_() -- wrong type of transform is provided (space2 != space4)");
  if (ipjq_tform->space1() != r12info_->refinfo()->docc_act())
    throw std::runtime_error("R12IntEval::compute_R_vbsneqobs_() -- wrong type of transform is provided (space1 != act_occ)");
  if (ipjq_tform->space2() != r12info_->refinfo()->orbs())
    throw std::runtime_error("R12IntEval::compute_R_vbsneqobs_() -- wrong type of transform is provided (space1 != all)");
  if (!ijpq_acc->is_committed())
    ipjq_tform->compute(intparams_);
  if (!ijpq_acc->is_active())
    ijpq_acc->activate();

  tim_enter("R intermediate");

  Ref<MessageGrp> msg = r12info_->msg();
  int me = msg->me();
  int nproc = msg->n();
  ExEnv::out0() << endl << indent
    << "Entered R intermediate evaluator" << endl;
  ExEnv::out0() << incindent;

  Ref<MOIndexSpace> act_occ_space = ipjq_tform->space1();
  Ref<MOIndexSpace> space2 = ipjq_tform->space2();
  Ref<MOIndexSpace> space4 = ipjq_tform->space4();
  const int rank2 = space2->rank();
  const int rank4 = space4->rank();

  MOPairIterFactory PIFactory;
  Ref<SpatialMOPairIter> ij_iter = PIFactory.mopairiter(act_occ_space,act_occ_space);
  Ref<SpatialMOPairIter> pq_iter = PIFactory.mopairiter(space2,space4);

  int nij_aa = ij_iter->nij_aa();          // Number of alpha-alpha pairs (i > j)
  int nij_ab = ij_iter->nij_ab();          // Number of alpha-beta pairs
  if (debug_) {
    ExEnv::out0() << indent << "nij_aa = " << nij_aa << endl;
    ExEnv::out0() << indent << "nij_ab = " << nij_ab << endl;
  }

  // Compute the number of tasks that have full access to the integrals
  // and split the work among them
  vector<int> proc_with_ints;
  int nproc_with_ints = tasks_with_ints_(ijpq_acc,proc_with_ints);

  for(ij_iter->start();int(*ij_iter.pointer());ij_iter->next()) {

    const int ij = ij_iter->ij();
    // Figure out if this task will handle this ij
    int ij_proc = ij%nproc_with_ints;
    if (ij_proc != proc_with_ints[me])
      continue;
    const int i = ij_iter->i();
    const int j = ij_iter->j();
    const int ij_aa = ij_iter->ij_aa();
    const int ij_ab = ij_iter->ij_ab();
    const int ji_ab = ij_iter->ij_ba();

    if (debug_)
      ExEnv::outn() << indent << "task " << me << ": working on (i,j) = " << i << "," << j << " " << endl;

    // Get (|1/r12|) integrals
    tim_enter("MO ints retrieve");
    double *ijxy_buf_r12 = ijpq_acc->retrieve_pair_block(i,j,corrfactor_->tbint_type_f12());
    double *jixy_buf_r12 = ijpq_acc->retrieve_pair_block(j,i,corrfactor_->tbint_type_f12());
    tim_exit("MO ints retrieve");

    if (debug_)
      ExEnv::outn() << indent << "task " << me << ": obtained ij blocks" << endl;

    for(pq_iter->start();int(*pq_iter.pointer());pq_iter->next()) {

      const int p = pq_iter->i();
      const int q = pq_iter->j();
      const int pq_aa = pq_iter->ij_aa();
      const int pq_ab = pq_iter->ij_ab();
      const int pq_ba = pq_iter->ij_ba();

      const int pq_offset = p*rank4 + q;
      const double r12_ipjq = ijxy_buf_r12[pq_offset];
      const double r12_jpiq = jixy_buf_r12[pq_offset];
      Rab.set_element(ij_ab,pq_ab,r12_ipjq);
      Rab.set_element(ji_ab,pq_ba,r12_ipjq);
      Rab.set_element(ij_ab,pq_ba,r12_jpiq);
      Rab.set_element(ji_ab,pq_ab,r12_jpiq);

      if (ij_aa != -1 && pq_aa != -1) {
        const double R_aa_ijpq = (r12_ipjq - r12_jpiq);
        Raa.set_element(ij_aa,pq_aa,R_aa_ijpq);
      }

    }

    ijpq_acc->release_pair_block(i,j,corrfactor_->tbint_type_f12());
    ijpq_acc->release_pair_block(j,i,corrfactor_->tbint_type_f12());
  }

  globally_sum_scmatrix_(Raa);
  globally_sum_scmatrix_(Rab);

  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Exited R intermediate evaluator" << endl;
  tim_exit("R intermediate");
}

void
R12IntEval::compute_amps_()
{
  if (Amps_.nonnull())
    return;

  Ref<MOIndexSpace> act_vir_space = r12info_->vir_act();
  Ref<MOIndexSpace> occ_space = r12info_->refinfo()->docc();
  Ref<MOIndexSpace> ribs_space = r12info_->ribs_space();

  MOPairIterFactory PIFactory;
  RefSCDimension ij_aa = dim_ij_aa_;
  RefSCDimension ij_ab = dim_ij_ab_;
  RefSCDimension vv_aa = PIFactory.scdim_aa(act_vir_space,act_vir_space);
  RefSCDimension vv_ab = PIFactory.scdim_ab(act_vir_space,act_vir_space);
  RefSCDimension oo_aa = PIFactory.scdim_aa(occ_space,occ_space);
  RefSCDimension oo_ab = PIFactory.scdim_ab(occ_space,occ_space);
  RefSCDimension ov_aa = PIFactory.scdim_aa(occ_space,act_vir_space);
  RefSCDimension ov_ab = PIFactory.scdim_ab(occ_space,act_vir_space);
  RefSCDimension ox_aa = PIFactory.scdim_aa(occ_space,ribs_space);
  RefSCDimension ox_ab = PIFactory.scdim_ab(occ_space,ribs_space);
  Ref<SCMatrixKit> kit = new LocalSCMatrixKit();

  RefSCMatrix T2_aa  = kit->matrix(ij_aa,vv_aa);
  RefSCMatrix T2_ab  = kit->matrix(ij_ab,vv_ab);
  RefSCMatrix Rvv_aa = kit->matrix(ij_aa,vv_aa);
  RefSCMatrix Rvv_ab = kit->matrix(ij_ab,vv_ab);
  RefSCMatrix Roo_aa = kit->matrix(ij_aa,oo_aa);
  RefSCMatrix Roo_ab = kit->matrix(ij_ab,oo_ab);
  RefSCMatrix Rvo_aa = kit->matrix(ij_aa,ov_aa);
  RefSCMatrix Rvo_ab = kit->matrix(ij_ab,ov_ab);
  RefSCMatrix Rxo_aa = kit->matrix(ij_aa,ox_aa);
  RefSCMatrix Rxo_ab = kit->matrix(ij_ab,ox_ab);
  compute_T2_vbsneqobs_();
  compute_R_vbsneqobs_(get_tform_("(ia|jb)"),Rvv_aa,Rvv_ab);
  compute_R_vbsneqobs_(get_tform_("(im|jn)"),Roo_aa,Roo_ab);
  compute_R_vbsneqobs_(get_tform_("(im|ja)"),Rvo_aa,Rvo_ab);
  compute_R_vbsneqobs_(get_tform_("(im|jy)"),Rxo_aa,Rxo_ab);

  Amps_ = new R12Amplitudes(T2_aa, T2_ab, Rvv_aa, Rvv_ab, Roo_aa, Roo_ab, Rvo_aa, Rvo_ab, Rxo_aa, Rxo_ab);
}

Ref<R12Amplitudes>
R12IntEval::amps()
{
  if (Amps_.null())
    compute_amps_();
  return Amps_;
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
