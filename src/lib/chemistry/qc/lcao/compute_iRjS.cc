//
// compute_iRjS.cc
//
// Copyright (C) 2008 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
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

#include<chemistry/qc/lcao/transform_iRjS.h>
#include<chemistry/qc/lcao/transform_13inds.h>
#include <util/misc/print.h>
#include <util/state/state.h>
#include <util/state/state_text.h>
#include <util/state/state_bin.h>
#include <math/scmat/local.h>

using namespace sc;

using namespace std;
using namespace sc;

#define SINGLE_THREAD_E13   0
#define PRINT2Q 0
#define PRINT_NUM_TE_TYPES 1

/*-------------------------------------
  Based on MBPT2::compute_mp2_energy()
 -------------------------------------*/
void
TwoBodyMOIntsTransform_iRjS::compute()
{
  const int rank1 = space1()->rank();
  const int rank2 = space2()->rank();
  const int rank3 = space3()->rank();
  const int rank4 = space4()->rank();

  init_acc();
  // if all integrals are already available -- do nothing
  if (restart_orbital_ == rank1)
    return;

  Ref<Integral> integral = factory_->integral();
  Ref<GaussianBasisSet> bs1 = space1()->basis();
  Ref<GaussianBasisSet> bs2 = space2()->basis();
  Ref<GaussianBasisSet> bs3 = space3()->basis();
  Ref<GaussianBasisSet> bs4 = space4()->basis();
  int nbasis1 = bs1->nbasis();
  int nbasis2 = bs2->nbasis();
  int nbasis3 = bs3->nbasis();
  int nbasis4 = bs4->nbasis();
  int nshell1 = bs1->nshell();
  int nshell2 = bs2->nshell();
  int nshell3 = bs3->nshell();
  int nshell4 = bs4->nshell();
  int nfuncmax1 = bs1->max_nfunction_in_shell();
  int nfuncmax2 = bs2->max_nfunction_in_shell();
  int nfuncmax3 = bs3->max_nfunction_in_shell();
  int nfuncmax4 = bs4->max_nfunction_in_shell();
  const size_t memgrp_blocksize = memgrp_blksize();

  // log2 of the erep tolerance
  // (erep < 2^tol => discard)
  const int tol = (int) (-10.0/log10(2.0));  // discard ints smaller than 10^-20

  int aoint_computed = 0;

  std::string tim_label("tbint_tform_iRjS ");
  tim_label += name_;
  Timer tim(tim_label);

  print_header();

  int me = msg_->me();
  int nproc = msg_->n();
  const int restart_orb = restart_orbital();
  const int nijmax = compute_nij(batchsize_,rank3,nproc,me);
  // allocate memory for MemoryGrp
  const size_t localmem = nijmax * num_te_types() * memgrp_blksize();
  alloc_mem(localmem);

  /////////////////////////////////////
  //  Begin transformation loops
  /////////////////////////////////////

  // debug print
  if (debug() >= DefaultPrintThresholds::fine) {
    ExEnv::outn() << indent
          << scprintf("node %i, begin loop over i-batches",me) << endl;
  }
  // end of debug print

  // Initialize the integrals
  integral->set_storage(max_memory_ - localmem);
  integral->set_basis(bs1,bs2,bs3,bs4);
  Ref<TwoBodyInt>* tbints = new Ref<TwoBodyInt>[thr_->nthread()];
  for (int i=0; i<thr_->nthread(); i++) {
    tbints[i] = tbintdescr_->inteval();
  }
  if (debug() >= DefaultPrintThresholds::diagnostics)
    ExEnv::out0() << indent << scprintf("Memory used for integral storage:       %i Bytes",
      integral->storage_used()) << endl;

  Ref<ThreadLock> lock = thr_->new_lock();
  TwoBodyMOIntsTransform_13Inds** e13thread = new TwoBodyMOIntsTransform_13Inds*[thr_->nthread()];
  for (int i=0; i<thr_->nthread(); i++) {
    e13thread[i] = new TwoBodyMOIntsTransform_13Inds(this,i,thr_->nthread(),lock,tbints[i],
                                                     this->log2_epsilon(),debug());
  }

  /*-----------------------------------
    Start the integrals transformation
   -----------------------------------*/
  Timer tim_pass("mp2-r12/a passes");
  if (me == 0 && top_mole_.nonnull() && top_mole_->if_to_checkpoint()) {
    StateOutBin stateout(top_mole_->checkpoint_file());
    SavableState::save_state(top_mole_,stateout);
    ExEnv::out0() << indent << "Checkpointed the wave function" << endl;
  }

  for (int pass=0; pass<npass_; pass++) {

    ExEnv::out0() << indent << "Beginning pass " << pass+1 << endl;

    int ni = batchsize_;
    int i_offset = restart_orb + pass*ni;
    if (pass == npass_ - 1)
      ni = rank1 - batchsize_*(npass_-1);

    // Compute number of of i,j pairs on each node during current pass for
    // two-el integrals
    int nij = compute_nij(ni,rank3,nproc,me);

    // debug print
    if (debug() >= DefaultPrintThresholds::fine)
      ExEnv::outn() << indent << "node " << me << ", nij = " << nij << endl;
    // end of debug print

    // Allocate and initialize some arrays
    // (done here to avoid having these arrays
    // overlap with arrays allocated later)

    // Allocate (and initialize) some arrays

    double* integral_ijRS = (double*) mem_->localdata();
    memset(integral_ijRS, 0, num_te_types()*nij*memgrp_blocksize);
    integral_ijRS = 0;
    mem_->sync();
    ExEnv::out0() << indent
          << scprintf("Begin loop over shells (ints, 1+2 q.t.)") << endl;

    // Do the two electron integrals and the first two quarter transformations
    Timer tim12("ints+1qt+2qt");
    shell_pair_data()->init();
    for (int i=0; i<thr_->nthread(); i++) {
      e13thread[i]->set_i_offset(i_offset);
      e13thread[i]->set_ni(ni);
      thr_->add_thread(i,e13thread[i]);
#     if SINGLE_THREAD_E13
      e13thread[i]->run();
#     endif
    }
#   if !SINGLE_THREAD_E13
    thr_->start_threads();
    thr_->wait_threads();
#   endif
    tim12.exit();
    ExEnv::out0() << indent << "End of loop over shells" << endl;

    mem_->sync();  // Make sure ijRS is complete on each node before continuing
    integral_ijRS = (double*) mem_->localdata();

#if PRINT2Q
    {
      for(int te_type=0; te_type<PRINT_NUM_TE_TYPES; te_type++) {
        for (int i = 0; i<ni; i++) {
          for (int q = 0; q<nbasis2; q++) {
            for (int j = 0; j<rank3; j++) {
              int ij = i*rank3+j;
              int ij_local = ij/nproc;
              if (ij%nproc == me) {
                const double* ijRS_ints = (const double*) ((size_t)integral_ijRS + (ij_local*num_te_types()+te_type)*memgrp_blocksize);
                for (int s = 0; s<nbasis4; s++) {
                  double value = ijRS_ints[s*rank2+q];
                  printf("2Q: type = %d (%d %d|%d %d) = %12.8f\n",
                         te_type,i+i_offset,q,j,s,value);
                }
              }
            }
          }
        }
      }
    }
#endif

    // Sync up tasks before integrals are committed
    mem_->sync();

    // Push locally stored integrals to an accumulator
    // This could involve storing the data to disk or simply remembering the pointer
    Timer tim_mostore("MO ints store");
    ints_acc_->activate();
    detail::store_memorygrp(ints_acc_,mem_,restart_orbital_,ni,memgrp_blocksize);
    if (ints_acc_->data_persistent()) ints_acc_->deactivate();
    // if didn't throw can safely update the counter
    restart_orbital_ += ni;
    tim_mostore.exit();
    mem_->sync();

    if (me == 0 && top_mole_.nonnull() && top_mole_->if_to_checkpoint()) {
      StateOutBin stateout(top_mole_->checkpoint_file());
      SavableState::save_state(top_mole_,stateout);
      ExEnv::out0() << indent << "Checkpointed the wave function" << endl;
    }

  } // end of loop over passes
  tim_pass.exit();

  for (int i=0; i<thr_->nthread(); i++) {
    delete e13thread[i];
  }
  delete[] e13thread;
  delete[] tbints; tbints = 0;

  tim.exit();

  if (me == 0 && top_mole_.nonnull() && top_mole_->if_to_checkpoint()) {
    StateOutBin stateout(top_mole_->checkpoint_file());
    SavableState::save_state(top_mole_,stateout);
    ExEnv::out0() << indent << "Checkpointed the wave function" << endl;
  }

#define DEBUG 0
#if DEBUG
        {
          ExEnv::out0() << indent << "TwoBodyMOIntsTransform_iRjS: name = " << this->name() << endl;
          ints_acc_->activate();
          for(int te_type=0; te_type<num_te_types(); ++te_type) {
            for(int i=0; i<rank1; ++i) {
              for(int j=0; j<rank3; ++j) {
                ExEnv::out0() << indent << "te_type = " << te_type << " i = " << i << " j = " << j << endl;
                Ref<SCMatrixKit> kit = new LocalSCMatrixKit;
                RefSCMatrix yx_buf_mat = kit->matrix(new SCDimension(rank4),
                                                     new SCDimension(rank2));

                const double* yx_buf = ints_acc_->retrieve_pair_block(i, j, te_type);
                yx_buf_mat.assign(yx_buf);
                yx_buf_mat.print("block");
                ints_acc_->release_pair_block(i, j, te_type);

              }
            }
          }
          if (ints_acc_->data_persistent()) ints_acc_->deactivate();
        }
#endif

  print_footer();

  // memory used by MemoryGrp can now be purged unless ints_acc_ uses it
  if (ints_acc_->data_persistent()) dealloc_mem();

}



/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
