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
#include <chemistry/qc/mbptr12/utils.h>

using namespace std;
using namespace sc;

void
R12IntEval::init_intermeds_g12_(SpinCase2 spincase)
{
  Ref<MessageGrp> msg = r12info()->msg();
  const int me = msg->me();

  tim_enter("\"diagonal\" part of G12 intermediates");
  ExEnv::out0() << endl << indent
	       << "Entered G12 diagonal intermediates evaluator" << endl;
  ExEnv::out0() << incindent;

  const unsigned int num_f12 = corrfactor()->nfunctions();
  const SpinCase1 spin1 = case1(spincase);
  const SpinCase1 spin2 = case2(spincase);
  /// compared to the needed (ik|jl), there's a chance (im|jn) has been done before (VBS != OBS?)
  const Ref<MOIndexSpace>& ispace = r12info()->refinfo()->occ_act(spin1);
  const Ref<MOIndexSpace>& mspace = r12info()->refinfo()->occ(spin1);
  const Ref<MOIndexSpace>& jspace = r12info()->refinfo()->occ_act(spin2);
  const Ref<MOIndexSpace>& nspace = r12info()->refinfo()->occ(spin2);
  Ref<MOIntsTransformFactory> tfactory = r12info()->tfactory();
  tfactory->set_spaces(ispace,mspace,jspace,nspace);
  std::vector<std::string> imjn_name;  imjn_name.resize(num_f12);
  std::vector<std::vector<std::string> > im2jn_name;  im2jn_name.resize(num_f12);
  // Compute transforms
  for(int f=0; f<num_f12; f++) {
    const std::string imjn_label = transform_label(ispace,mspace,jspace,nspace,f);
    imjn_name[f] = imjn_label;

    Ref<TwoBodyMOIntsTransform> imjn_tform = tform_map_[imjn_label];
    if (imjn_tform.null()) {
      imjn_tform = tfactory->twobody_transform_13(imjn_label,corrfactor()->callback());
      imjn_tform->set_num_te_types(corrfactor()->num_tbint_types());
      tform_map_[imjn_label] = imjn_tform;
      Ref<IntParams> params = new IntParamsG12(corrfactor()->function(f),LinearR12::CorrelationFactor::zero_exponent_geminal());
      imjn_tform->compute(params);
      }

    // second loop over correlation functions
    //
    // 2) get (im|jn) integrals of [g12,[t1,g12]] and g12*g12 operator (use integrals with the exponent multiplied by 2, and additionally [g12,[t1,g12]] integral needs to be scaled by 0.25 to take into account that real exponent is half what the integral library thinks)
    //    these integrals used to compute X and B
    for(int g=0; g<=f; g++) {
      // (im|2|jn) name
      const std::string im2jn_label = transform_label(ispace,mspace,jspace,nspace,f,g);;
      im2jn_name[f].push_back(im2jn_label);
      Ref<TwoBodyMOIntsTransform> im2jn_tform = tform_map_[im2jn_label];
      if (im2jn_tform.null()) {
        im2jn_tform = tfactory->twobody_transform_13(im2jn_label,corrfactor()->callback());
        im2jn_tform->set_num_te_types(corrfactor()->num_tbint_types());
        tform_map_[im2jn_label] = im2jn_tform;
        Ref<IntParams> params = new IntParamsG12(corrfactor()->function(f),
                                                 corrfactor()->function(g));
        im2jn_tform->compute(params);
      }
    }
  }

  const int ni = ispace->rank();
  const int nj = jspace->rank();
  const int nm = mspace->rank();
  const int nn = nspace->rank();
  const int nij = ni*nj;
  const int nfzc = r12info()->refinfo()->nfzc();

    // If same spin -- will compute different spin and antisymmetrize at the end
  const bool need_to_antisymmetrize = (spincase != AlphaBeta && ispace == jspace);
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

  
  ExEnv::out0() << indent << "Begin computation of intermediates" << endl;
  tim_enter("intermediates");
  SpatialMOPairIter_neq ij_iter(ispace,jspace);
  SpatialMOPairIter_neq kl_iter(ispace,jspace);
  const int k_offset = nm - ni;
  const int l_offset = nn - nj;
  if (nfzc != k_offset || nfzc != l_offset)
    throw ProgrammingError("R12IntEval::init_intermeds_g12_() -- given orbital spaces do not match",__FILE__,__LINE__);
  
  // Compute the number of tasks that have full access to the integrals
  // and split the work among them
  // WARNING Assuming here that all accumulators have the same availability, otherwise
  // things will get messy
  vector<int> proc_with_ints;
  Ref<R12IntsAcc> any_accum = get_tform_(imjn_name[0])->ints_acc();
  int nproc_with_ints = tasks_with_ints_(any_accum,proc_with_ints);

  //////////////////////////////////////////////////////////////
  //
  // Evaluation of the intermediates proceeds as follows:
  //
  // loop over contracted geminals f
  //   loop over ij
  //     load (ij|f12[f]/r12|kl) into memory
  //     compute V[ij][kl]
  //
  //     loop over contracted geminals g
  //       load (ij|f12[f]*f12[g]|kl), (ij| [f12[f],[T1,f12[g]]] |kl)
  //       compute X[ij][kl] and B[ij][kl]
  //     end g loop
  //
  //   end ij loop
  // end f loop
  //
  /////////////////////////////////////////////////////////////////////////////////

  if (any_accum->has_access(me)) {

    // f loop
    for(int f=0; f<num_f12; f++) {
      const int f_offset = f*nij;
      Ref<R12IntsAcc> ijmn_acc = get_tform_(imjn_name[f])->ints_acc();
      
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
        const int ij_abs = ij_ab + f_offset;
        
        if (debug_)
          ExEnv::outn() << indent << "task " << me << ": working on (i,j) = " << i << "," << j << " " << endl;

        //
        // Compute V
        //
        tim_enter("MO ints retrieve");
        const double *ijmn_buf_f12eri = ijmn_acc->retrieve_pair_block(i,j,corrfactor()->tbint_type_f12eri());
        tim_exit("MO ints retrieve");
        if (debug_)
          ExEnv::outn() << indent << "task " << me << ": obtained ij blocks for V" << endl;
        // kl loop
        for(kl_iter.start();int(kl_iter);kl_iter.next()) {
          const int k = kl_iter.i();
          const int kk = k + k_offset;
          const int l = kl_iter.j();
          const int ll = l + l_offset;
          const int kl_ab = kl_iter.ij_ab();
          const int kkll = kk*nn + ll;

          const double V_ijkl = ijmn_buf_f12eri[kkll];
          V.accumulate_element(ij_abs,kl_ab,V_ijkl);
        }
        ijmn_acc->release_pair_block(i,j,corrfactor()->tbint_type_f12eri());

        //
        // Compute X and B
        //
        
        // g loop
        for(int g=0; g<=f; g++) {
          const int g_offset = g*nij;
          const bool f_eq_g = (f == g);
          Ref<R12IntsAcc> ij2mn_acc = get_tform_(im2jn_name[f][g])->ints_acc();
          
          tim_enter("MO ints retrieve");
          const double *ijmn_buf_f12f12 = ij2mn_acc->retrieve_pair_block(i,j,corrfactor()->tbint_type_f12f12());
          const double *ijmn_buf_f12t1f12 = ij2mn_acc->retrieve_pair_block(i,j,corrfactor()->tbint_type_f12t1f12());
          tim_exit("MO ints retrieve");
          if (debug_)
            ExEnv::outn() << indent << "task " << me << ": obtained ij blocks for X and B" << endl;

          // kl loop
          for(kl_iter.start();int(kl_iter);kl_iter.next()) {
            const int k = kl_iter.i();
            const int kk = k + k_offset;
            const int l = kl_iter.j();
            const int ll = l + l_offset;
            const int kl_ab = kl_iter.ij_ab();
            const int kkll = kk*nn + ll;
            const int kl_abs = kl_ab + g_offset;

            const double X_ijkl = ijmn_buf_f12f12[kkll];
            const double B_ijkl = ijmn_buf_f12t1f12[kkll];
            X.accumulate_element(ij_abs,kl_abs,X_ijkl);
            B.accumulate_element(ij_abs,kl_abs,B_ijkl);
            if (!f_eq_g) {
              X.accumulate_element(kl_abs,ij_abs,X_ijkl);
              B.accumulate_element(kl_abs,ij_abs,B_ijkl);
            }
          } // end of kl loop
          ij2mn_acc->release_pair_block(i,j,corrfactor()->tbint_type_f12f12());
          ij2mn_acc->release_pair_block(i,j,corrfactor()->tbint_type_f12t1f12());
        } // end of g loop
      } // end of ij loop
    } // end of f loop
  } // if (have_access) block

  // Tasks that don't do any work here still need to create these timers
  tim_enter("MO ints retrieve");
  tim_exit("MO ints retrieve");

  tim_exit("intermediates");
  ExEnv::out0() << indent << "End of computation of intermediates" << endl;
  for(int f=0; f<num_f12; f++) {
    Ref<R12IntsAcc> ijmn_acc = get_tform_(imjn_name[f])->ints_acc();
    ijmn_acc->deactivate();
    for(int g=0; g<=f; g++) {
      Ref<R12IntsAcc> ij2mn_acc = get_tform_(im2jn_name[f][g])->ints_acc();
      ij2mn_acc->deactivate();
    }
  }

  if (need_to_antisymmetrize) {
    // antisymmetrize I and add to I_[spincase]
    antisymmetrize(V_[spincase],V,ispace,ispace,true);
    antisymmetrize(X_[spincase],X,ispace,ispace,true);
    antisymmetrize(B_[spincase],B,ispace,ispace,true);
  }
  V = 0; X = 0; B = 0;
  
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
