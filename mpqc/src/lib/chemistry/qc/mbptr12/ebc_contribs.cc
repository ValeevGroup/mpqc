//
// ebc_contribs.cc
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
#include <chemistry/qc/mbptr12/print_scmat_norms.h>
#include <chemistry/qc/mbptr12/creator.h>
#include <chemistry/qc/mbptr12/container.h>
#include <chemistry/qc/mbptr12/utils.h>
#include <chemistry/qc/mbptr12/utils.impl.h>

using namespace std;
using namespace sc;

#define TEST_T2 0
#define TEST_A 0
// if set to 1 then use f+k rather than f to compute A
#define A_DIRECT_EXCLUDE_K 0

//
// these are for testing purposes only
//
// use the commutator form A
#define USE_A_COMM_IN_B_EBC 0
#define ACOMM_INCLUDE_TR_ONLY 0
#define ACOMM_INCLUDE_R_ONLY 0


void
R12IntEval::compute_A_direct_(RefSCMatrix& A,
                              const Ref<MOIndexSpace>& space1,
                              const Ref<MOIndexSpace>& space2,
                              const Ref<MOIndexSpace>& space3,
                              const Ref<MOIndexSpace>& space4,
                              const Ref<MOIndexSpace>& fspace2,
                              const Ref<MOIndexSpace>& fspace4)
{
  // are particles 1 and 2 equivalent?
  const bool part1_equiv_part2 = (space1==space3 && space2 == space4);
  
  const unsigned int nf12 = corrfactor()->nfunctions();
  // create transforms, if needed
  std::vector< Ref<TwoBodyIntDescr> > descrs; // get 1 3 |F12| 2 4_f
  TwoBodyIntDescrCreator descr_creator(corrfactor(),
                                       r12info()->integral(),
                                       true,false);
  fill_container(descr_creator,descrs);
  
  tim_enter("A intermediate (direct)");
  std::ostringstream oss;
  oss << "<" << space1->id() << " " << space3->id() << "|A|"
      << space2->id() << " " << space4->id() << ">";
  const std::string label = oss.str();
  ExEnv::out0() << endl << indent
                << "Entered \"direct\" A intermediate (" << label << ") evaluator" << endl
                << incindent;
  //
  // ij|A|kl = ij|f12|kl_f, symmetrized if part_equiv_part2
  //
  std::vector< Ref<TwoBodyMOIntsTransform> > tforms4f; // get 1 3 |F12| 2 4_f
  compute_F12_(A,space1,space2,space3,fspace4,tforms4f,descrs);
  if (part1_equiv_part2) {
    symmetrize<true,false,false>(A,A,space1,space2);
  }
  else {
    std::vector< Ref<TwoBodyMOIntsTransform> > tforms2f;
    compute_F12_(A,space1,fspace2,space3,space4,tforms2f,descrs);
    A.scale(0.5);
  }

  ExEnv::out0() << decindent << indent << "Exited \"direct\" A intermediate (" << label << ") evaluator" << endl;
  tim_exit("A intermediate (direct)");
}


#if 0
void
R12IntEval::compute_A_via_commutator_()
{
  if (evaluated_)
    return;

  // This functions assumes that virtuals are expanded in the same basis
  // as the occupied orbitals
  if (!r12info_->basis_vir()->equiv(r12info_->basis()))
    throw std::runtime_error("R12IntEval::compute_A_via_commutator_() -- should not be called when the basis set for virtuals \
differs from the basis set for occupieds");

  Ref<TwoBodyMOIntsTransform> ipjq_tform = get_tform_("(ip|jq)");
  Ref<R12IntsAcc> ijpq_acc = ipjq_tform->ints_acc();
  if (!ijpq_acc->is_committed())
    ipjq_tform->compute();
  if (!ijpq_acc->is_active())
    ijpq_acc->activate();

  tim_enter("A intermediate via [T,r]");

  Ref<MessageGrp> msg = r12info_->msg();
  int me = msg->me();
  int nproc = msg->n();
  ExEnv::out0() << endl << indent
    << "Entered A amplitude (via [T,r]) evaluator" << endl;
  ExEnv::out0() << incindent;

  const int noso = r12info_->mo_space()->rank();
  const int nocc = r12info_->nocc();
  const int nfzv = r12info_->nfzv();
  Ref<MOIndexSpace> mo_space = r12info_->obs_space();
  Ref<MOIndexSpace> act_occ_space = r12info_->act_occ_space();
  Ref<MOIndexSpace> act_vir_space = r12info_->act_vir_space();

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
  int nproc_with_ints = tasks_with_ints_(ijpq_acc,proc_with_ints);

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
    double *ijxy_buf_r12 = ijpq_acc->retrieve_pair_block(i,j,R12IntsAcc::r12);
    double *ijxy_buf_r12t1 = ijpq_acc->retrieve_pair_block(i,j,R12IntsAcc::r12t1);
    double *jixy_buf_r12t1 = ijpq_acc->retrieve_pair_block(j,i,R12IntsAcc::r12t1);
    tim_exit("MO ints retrieve");

    if (debug_)
      ExEnv::outn() << indent << "task " << me << ": obtained ij blocks" << endl;

    RefDiagSCMatrix act_occ_evals = r12info_->act_occ_space()->evals();
    RefDiagSCMatrix all_evals = r12info_->obs_space()->evals();
    double R_aa_ijab = 0.0;
    double R_ab_ijab = 0.0;
    double TR_aa_ijab = 0.0;
    double TR_ab_ijab = 0.0;

    for(ab_iter.start();int(ab_iter);ab_iter.next()) {

      const int a = ab_iter.i();
      const int b = ab_iter.j();
      const int ab_aa = ab_iter.ij_aa();
      const int ab_ab = ab_iter.ij_ab();
      const int ba_ab = ab_iter.ij_ba();

      const int aa = a + nocc;
      const int bb = b + nocc;

      const int ab_offset = aa*noso+bb;
      const int ba_offset = bb*noso+aa;
      const double r12_iajb = ijxy_buf_r12[ab_offset];
      const double r12_ibja = ijxy_buf_r12[ba_offset];
      const double t1r12_iajb = -ijxy_buf_r12t1[ab_offset];
      const double t1r12_ibja = -ijxy_buf_r12t1[ba_offset];
      const double t2r12_iajb = -jixy_buf_r12t1[ba_offset];
      const double t2r12_ibja = -jixy_buf_r12t1[ab_offset];
#if ACOMM_INCLUDE_TR_ONLY
      double Aab_ij_ab = 0.5 * ( -(t1r12_iajb + t2r12_iajb) );
      double Aab_ij_ba = 0.5 * ( -(t1r12_ibja + t2r12_ibja) );
#elif ACOMM_INCLUDE_R_ONLY
      double Aab_ij_ab = 0.5 * ( r12_iajb );
      double Aab_ij_ba = 0.5 * ( r12_ibja );
#else
      double Aab_ij_ab = 0.5 * ( -(t1r12_iajb + t2r12_iajb) - (all_evals(aa) + all_evals(bb) -
                                                              act_occ_evals(i) - act_occ_evals(j))*r12_iajb );
      double Aab_ij_ba = 0.5 * ( -(t1r12_ibja + t2r12_ibja) - (all_evals(aa) + all_evals(bb) -
                                                              act_occ_evals(i) - act_occ_evals(j))*r12_ibja );
#endif
      Ac_ab_.set_element(ij_ab,ab_ab,Aab_ij_ab);
      Ac_ab_.set_element(ji_ab,ba_ab,Aab_ij_ab);
      Ac_ab_.set_element(ji_ab,ab_ab,Aab_ij_ba);
      Ac_ab_.set_element(ij_ab,ba_ab,Aab_ij_ba);

      if (ij_aa != -1 && ab_aa != -1) {
        Ac_aa_.set_element(ij_aa,ab_aa,Aab_ij_ab - Aab_ij_ba);
      }

    }

    ijpq_acc->release_pair_block(i,j,R12IntsAcc::r12);
    ijpq_acc->release_pair_block(i,j,R12IntsAcc::r12t1);
    ijpq_acc->release_pair_block(j,i,R12IntsAcc::r12t1);
  }

  globally_sum_intermeds_();

#if TEST_A
  Ac_aa_.print("Alpha-alpha A intermediate");
  Ac_ab_.print("Alpha-beta A intermediate");
#endif

  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Exited A amplitude (via [T,r]) evaluator" << endl;

  tim_exit("A intermediate via [T,r]");
}
#endif // Commented-out old code

void
R12IntEval::AT2_contrib_to_V_()
{
  if (evaluated_)
    return;
  if (r12info_->msg()->me() == 0) {
    for(unsigned int s=0; s<nspincases2(); s++) {
      SpinCase2 spin = static_cast<SpinCase2>(s);
      RefSCMatrix V = A_[s]*T2_[s].t();
      if (debug_ > 0) {
        std::string label = prepend_spincase(spin,"AT2 contribution to V");
        print_scmat_norms(V,label.c_str());
      }
    }
  }
  globally_sum_intermeds_();
}

void
R12IntEval::AF12_contrib_to_B_()
{
  if (evaluated_)
    return;
  if (r12info_->msg()->me() == 0) {
    for(unsigned int s=0; s<nspincases2(); s++) {
      SpinCase2 spin = static_cast<SpinCase2>(s);
      RefSCMatrix AR = A_[s]*F12_[s].t();
      const double scale = -0.5;
#if 0
#if USE_A_COMM_IN_B_EBC
      RefSCMatrix AR_aa = Ac_aa_*Raa_.t();
      RefSCMatrix AR_ab = Ac_ab_*Rab_.t();
      double scale = -0.5;
#else
      RefSCMatrix AR_aa = Aaa_*Raa_.t();
      RefSCMatrix AR_ab = Aab_*Rab_.t();
      double scale = -0.5;
#endif
#endif
      RefSCMatrix B = B_[s].clone();  B.assign(0.0);
      AR.scale(scale); B.accumulate(AR);
      RefSCMatrix ARt = AR.t();
      B.accumulate(ARt);

      const std::string label = prepend_spincase(spin,"B^{EBC} contribution");
      if (debug_ > 1) {
        B.print(label.c_str());
      }
      B_[s].accumulate(B);
      if (debug_ > 0) {
        print_scmat_norms(B,label.c_str());
      }
    }
  }
  globally_sum_intermeds_();
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
