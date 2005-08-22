//
// compute_vxb_a_new.cc
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

#include <sstream>
#include <stdlib.h>
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
#include <chemistry/qc/mbptr12/twoparticlecontraction.h>
#include <chemistry/qc/mbptr12/utils.h>

using namespace std;
using namespace sc;

#define PRINT_R12_INTERMED 0

using LinearR12::TwoParticleContraction;

namespace {
void permsymm_block(const RefSCMatrix& A, const Ref<MOIndexSpace>& ispace,
                    unsigned int row_offset, unsigned int col_offset);
}

/**
   R12IntEval::contrib_to_VXB_a_new_() computes V, X, and B contributions
   from <ij|xy> integrals
*/
void
R12IntEval::contrib_to_VXB_a_new_(const Ref<MOIndexSpace>& ispace,
                                  const Ref<MOIndexSpace>& xspace,
                                  const Ref<MOIndexSpace>& jspace,
                                  const Ref<MOIndexSpace>& yspace,
                                  SpinCase2 spincase,
                                  const Ref<TwoParticleContraction>& tpcontract)
{
  if (evaluated_)
    return;
  LinearR12::ABSMethod abs_method = r12info_->abs_method();

  tim_enter("mp2-f12a intermeds (new)");

  const unsigned int num_f12 = corrfactor()->nfunctions();
  Ref<MOIntsTransformFactory> tfactory = r12info()->tfactory();
  tfactory->set_spaces(ispace,xspace,jspace,yspace);
  std::vector<std::string> ixjy_name;
  for(int f12=0; f12<num_f12; f12++) {
    const std::string tform_name = transform_label(ispace,xspace,jspace,yspace,f12);
    ixjy_name.push_back(tform_name);
    Ref<TwoBodyMOIntsTransform> ixjy_tform = tform_map_[tform_name];
    if (ixjy_tform.null()) {
      ixjy_tform = tfactory->twobody_transform_13(tform_name,corrfactor()->callback());
      ixjy_tform->set_num_te_types(corrfactor()->num_tbint_types());
      tform_map_[tform_name] = ixjy_tform;
      Ref<IntParams> params = new IntParamsG12(corrfactor()->function(f12),LinearR12::CorrelationFactor::zero_exponent_geminal());
      ixjy_tform->compute(params);
    }
  }
  const int ni = ispace->rank();
  const int nj = jspace->rank();
  const int nx = xspace->rank();
  const int ny = yspace->rank();
  const int nxy = nx*ny;
  const int nij = ni*nj;

  // If same spin -- will compute different spin and antisymmetrize at the end
  const bool need_to_antisymmetrize = (spincase != AlphaBeta && ispace == jspace);
  // still, need to check dimensions
  if (need_to_antisymmetrize) {
    const int nij_target = V_[spincase].coldim().n();
    if (nij_target != ni*(ni-1)/2)
      throw ProgrammingError("R12IntEval::contrib_to_VXB_a_new_() -- dimension of spaces 1 and 3 don't match requested spincase",
                             __FILE__,__LINE__);
    const int nf12_target = V_[spincase].rowdim().n();
    if (nf12_target != num_f12 * ni*(ni-1)/2)
      throw ProgrammingError("R12IntEval::contrib_to_VXB_a_new_() -- dimension of spaces 1 and 3 don't match requested spincase",
                             __FILE__,__LINE__);
  }
  const bool need_to_symmetrize = (ispace == jspace && xspace != yspace);
  const double perm_pfac = (need_to_symmetrize ? 2.0 : 1.0);

  RefSCMatrix V, X, B;
  if (!need_to_antisymmetrize) {
    V = V_[spincase];
    X = X_[spincase];
    B = B_[spincase];
  }
  else {
    Ref<SCMatrixKit> kit = V_[spincase].kit();
    RefSCDimension f12ab_dim = new SCDimension(num_f12 * nij);
    RefSCDimension ooab_dim = new SCDimension(nij);
    V = kit->matrix(f12ab_dim,ooab_dim);
    X = kit->matrix(f12ab_dim,f12ab_dim);
    B = kit->matrix(f12ab_dim,f12ab_dim);
    V.assign(0.0);
    X.assign(0.0);
    B.assign(0.0);
  }
  
  ExEnv::out0() << endl << indent
                << "Computing contribution to MP2-F12/A (GEBC) intermediates from "
                << ixjy_name[0] << " integrals" << endl;
  ExEnv::out0() << incindent;

  // Carry out the AO->MO transform and collect all accumulators
  std::vector< Ref<R12IntsAcc> > ijxy_acc;
  for(int f12=0; f12<num_f12; f12++) {
    Ref<TwoBodyMOIntsTransform> ixjy_tform = get_tform_(ixjy_name.at(f12));
    Ref<R12IntsAcc> acc = ixjy_tform->ints_acc();
    if (acc.null() || !acc->is_committed()) {
      Ref<IntParams> params = new IntParamsG12(corrfactor()->function(f12),LinearR12::CorrelationFactor::zero_exponent_geminal());
      ixjy_tform->compute(params);
    }
    if (!acc->is_active())
      acc->activate();
    ijxy_acc.push_back(acc);
  }

  /*--------------------------------
    Compute MP2-R12/A intermediates
    and collect on node0
   --------------------------------*/
  ExEnv::out0() << indent << "Begin computation of intermediates" << endl;
  tim_enter("intermediates");
  SpatialMOPairIter_neq ij_iter(ispace,jspace);
  SpatialMOPairIter_neq kl_iter(ispace,jspace);

  // Compute intermediates
  if (debug_)
    ExEnv::out0() << indent << "Ready to compute MP2-F12/A (GEBC) intermediates" << endl;

  // Compute the number of tasks that have full access to the integrals
  // and split the work among them
  // WARNING Assuming here that all accumulators have the same availability, otherwise
  // things will get messy
  vector<int> proc_with_ints;
  const int nproc_with_ints = tasks_with_ints_(ijxy_acc.at(0),proc_with_ints);

  //////////////////////////////////////////////////////////////
  //
  // Evaluation of the intermediates proceeds as follows:
  //
  // loop over contracted geminals f
  //   loop over ij
  //     load (ij|f12[f]|xy) (aka ij-sets) into memory
  //
  //     loop over kl
  //       load (kl|1/r12|xy) into memory
  //       compute V[ij][kl]
  //
  //     loop over contracted geminals g
  //       loop over kl
  //         load (kl|f12[g]|xy), (kl| [T1,f12[g]] |xy), and (lk| [T2,f12[g]] |xy)
  //           into memory
  //         compute X[ij][kl] and B[ij][kl]
  //       end kl loop
  //     end g loop
  //
  //   end ij loop
  // end f loop
  //
  /////////////////////////////////////////////////////////////////////////////////

  int me = r12info()->msg()->me();
  if (ijxy_acc[0]->has_access(me)) {

    // f loop
    for(int f=0; f<num_f12; f++) {
      const int f_offset = f*nij;
      // ij loop
      for(ij_iter.start();int(ij_iter);ij_iter.next()) {
        
        const int ij = ij_iter.ij();
        // Figure out if this task will handle this ij
        int ij_proc = ij%nproc_with_ints;
        if (ij_proc != proc_with_ints[me])
          continue;
        const int i = ij_iter.i();
        const int j = ij_iter.j();
        const int ij_ab = ij_iter.ij_ab();
        
        if (debug_)
          ExEnv::outn() << indent << "task " << me << ": working on (i,j) = " << i << "," << j << " " << endl;

        tim_enter("MO ints retrieve");
        const double *ijxy_buf_f12 = ijxy_acc[f]->retrieve_pair_block(i,j,corrfactor()->tbint_type_f12());
        tim_exit("MO ints retrieve");
        if (debug_)
          ExEnv::outn() << indent << "task " << me << ": obtained ij blocks" << endl;

        //////////////
        // Compute V
        //////////////
        // kl loop
        for(kl_iter.start(ij+1);int(kl_iter);kl_iter.next()) {
          
          const int kl = kl_iter.ij();
          // Figure out if this task will handle this kl
          int kl_proc = kl%nproc_with_ints;
          if (kl_proc != proc_with_ints[me])
            continue;
          const int k = kl_iter.i();
          const int l = kl_iter.j();
          const int kl_ab = kl_iter.ij_ab();
          
          if (debug_)
            ExEnv::outn() << indent << "task " << me << ": working on (k,l) = " << k << "," << l << " " << endl;
          
          // Get (|1/r12|) integrals
          tim_enter("MO ints retrieve");
          const double *klxy_buf_eri = ijxy_acc[0]->retrieve_pair_block(k,l,corrfactor()->tbint_type_eri());
          tim_exit("MO ints retrieve");
          if (debug_)
            ExEnv::outn() << indent << "task " << me << ": obtained kl blocks" << endl;
          
          tim_enter("MO ints contraction");
          
          double V_ijkl = tpcontract->contract(ijxy_buf_f12,klxy_buf_eri);
          V_ijkl *= perm_pfac;
          V.accumulate_element(f_offset+ij_ab,kl_ab,V_ijkl);
          if (debug_)
            ExEnv::out0() << " f = " << f << " ij = " << ij_ab << " kl = " << kl_ab << " V = " << V_ijkl << endl;
          tim_exit("MO ints contraction");
          
          ijxy_acc[0]->release_pair_block(k,l,corrfactor()->tbint_type_eri());
        } // end of kl loop
        
        ////////////////////
        // Compute X and B
        ////////////////////
        // g loop
        for(int g=0; g<num_f12; g++) {
          const int g_offset = g*nij;
          // kl loop
          for(kl_iter.start(ij+1);int(kl_iter);kl_iter.next()) {
            
            const int kl = kl_iter.ij();
            // Figure out if this task will handle this kl
            int kl_proc = kl%nproc_with_ints;
            if (kl_proc != proc_with_ints[me])
              continue;
            const int k = kl_iter.i();
            const int l = kl_iter.j();
            const int kl_ab = kl_iter.ij_ab();
            
            if (debug_)
              ExEnv::outn() << indent << "task " << me << ": working on (k,l) = " << k << "," << l << " " << endl;
            
            // Get (|f12|) and (|[T_i,f12]|) integrals
            tim_enter("MO ints retrieve");
            const double *klxy_buf_f12 = ijxy_acc[g]->retrieve_pair_block(k,l,corrfactor()->tbint_type_f12());
            const double *klxy_buf_t1f12 = ijxy_acc[g]->retrieve_pair_block(k,l,corrfactor()->tbint_type_t1f12());
            const double *klxy_buf_t2f12 = ijxy_acc[g]->retrieve_pair_block(k,l,corrfactor()->tbint_type_t2f12());
            tim_exit("MO ints retrieve");
            if (debug_)
              ExEnv::outn() << indent << "task " << me << ": obtained kl blocks" << endl;
            
            tim_enter("MO ints contraction");
            double X_ijkl = tpcontract->contract(ijxy_buf_f12,klxy_buf_f12);
            X_ijkl *= perm_pfac;
            X.accumulate_element(f_offset+ij_ab,g_offset+kl_ab,X_ijkl);
            if (debug_)
              ExEnv::out0() << " f = " << f << " ij = " << ij_ab
                            << " g = " << g << " kl = " << kl_ab
                            << " X = " << X_ijkl << endl;
            double B_ijkl = tpcontract->contract(ijxy_buf_f12,klxy_buf_t1f12);
            B_ijkl += tpcontract->contract(ijxy_buf_f12,klxy_buf_t2f12);
            B_ijkl *= perm_pfac;
            B.accumulate_element(f_offset+ij_ab,g_offset+kl_ab,B_ijkl);
            tim_exit("MO ints contraction");
            
            ijxy_acc[g]->release_pair_block(k,l,corrfactor()->tbint_type_f12());
            ijxy_acc[g]->release_pair_block(k,l,corrfactor()->tbint_type_t1f12());
            ijxy_acc[g]->release_pair_block(k,l,corrfactor()->tbint_type_t2f12());
          }
        }
          
        ijxy_acc[f]->release_pair_block(i,j,corrfactor()->tbint_type_f12());
      } // end of ij loop
    } // end of f loop
  } // end of (if I have access to integrals) block

  // Tasks that don't do any work here still need to create these timers
  tim_enter("MO ints retrieve");
  tim_exit("MO ints retrieve");
  tim_enter("MO ints contraction");
  tim_exit("MO ints contraction");

  tim_exit("intermediates");
  ExEnv::out0() << indent << "End of computation of intermediates" << endl;
  for(int f12=0; f12<num_f12; f12++) {
    ijxy_acc[f12]->deactivate();
  }
  
  // Symmetrize X and B intermediates with respect to geminal permutations
  for(int f=1; f<num_f12; f++) {
    const int f_off = f*nij;
    for(int g=0; g<f; g++) {
      const int g_off = g*nij;
      for(int ij=0; ij<nij; ij++) {
        const int IJf = f_off + ij;
        const int IJg = g_off + ij;
        for(int kl=0; kl<nij; kl++) {
          const int KLf = f_off + kl;
          const int KLg = g_off + kl;
          const double x = 0.5 * (X.get_element(IJf,KLg) + X.get_element(IJg,KLf));
          X.set_element(IJf,KLg,x);
          X.set_element(IJg,KLf,x);
          const double b = 0.5 * (B.get_element(IJf,KLg) + B.get_element(IJg,KLf));
          B.set_element(IJf,KLg,b);
          B.set_element(IJg,KLf,b);
        }
      }
    }
  }

  // Symmetrize B intermediate with respect to bra-ket permutation
  const int f12dim = B.coldim().n();
  for(int ij=0;ij<f12dim;ij++)
    for(int kl=0;kl<=ij;kl++) {
      double belem = 0.5*(B.get_element(ij,kl) + B.get_element(kl,ij));
      B.set_element(ij,kl,belem);
      B.set_element(kl,ij,belem);
    }
  
  // symmetrize all intermediates with respect to permutation of particles, if needed
  if (need_to_symmetrize) {
    for(int f=0; f<num_f12; f++) {
      const int f_off = f*nij;
      permsymm_block(V,ispace,f_off,0);
      for(int g=0; g<num_f12; g++) {
        const int g_off = g*nij;
        permsymm_block(X,ispace,f_off,g_off);
        permsymm_block(B,ispace,f_off,g_off);
      }
    }
  }

  if (need_to_antisymmetrize) {
    // antisymmetrize I and add to I_[spincase]
    antisymmetrize(V_[spincase],V,ispace,ispace,true);
    antisymmetrize(X_[spincase],X,ispace,ispace,true);
    antisymmetrize(B_[spincase],B,ispace,ispace,true);
  }
  
  globally_sum_intermeds_();
  
  ExEnv::out0() << decindent;
  ExEnv::out0() << indent
                << "Finished computing contribution to MP2-F12/A (GEBC) intermediates from "
                << ixjy_name[0] << " integrals" << endl << endl;

  tim_exit("mp2-f12a intermeds (new)");
  checkpoint_();
}

namespace {
void permsymm_block(const RefSCMatrix& A, const Ref<MOIndexSpace>& ispace,
                    unsigned int row_offset, unsigned int col_offset)
{
  SpatialMOPairIter_eq ij_iter(ispace);
  SpatialMOPairIter_eq kl_iter(ispace);
  // ij loop
  for(ij_iter.start();int(ij_iter);ij_iter.next()) {
    const int ij = ij_iter.ij_ab() + row_offset;
    const int ji = ij_iter.ij_ba() + row_offset;
    // kl loop
    for(kl_iter.start();int(kl_iter);kl_iter.next()) {
      const int kl = kl_iter.ij_ab() + col_offset;
      const int lk = kl_iter.ij_ba() + col_offset;
      const double A_ijkl = A.get_element(ij,kl);
      const double A_jilk = A.get_element(ji,lk);
      const double A1 = 0.5*(A_ijkl + A_jilk);
      A.set_element(ij,kl,A1);
      A.set_element(ji,lk,A1);
      const double A_ijlk = A.get_element(ij,lk);
      const double A_jikl = A.get_element(ji,kl);
      const double A2 = 0.5*(A_ijlk + A_jikl);
      A.set_element(ij,lk,A2);
      A.set_element(ji,kl,A2);
    }
  }
}

}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
