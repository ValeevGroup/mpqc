//
// transform_ijxy.cc
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

#ifdef __GNUG__
#pragma implementation
#endif

#include <stdexcept>

#include <util/misc/formio.h>
#include <util/state/state_bin.h>
#include <util/ref/ref.h>
#include <math/scmat/local.h>
#include <chemistry/qc/mbptr12/transform_ijxy.h>

// set to 1 when finished rewriting R12IntsAcc_MemoryGrp
#define HAVE_R12IA_MEMGRP 0
#if HAVE_R12IA_MEMGRP
  #include <chemistry/qc/mbptr12/r12ia_memgrp.h>
#endif
#include <cassert>

#include <chemistry/qc/mbptr12/r12ia_node0file.h>

// set to 1 when finished rewriting R12IntsAcc_MPIIO
#define HAVE_R12IA_MPIIO 0
#ifdef HAVE_MPIIO
#  if HAVE_R12IA_MPIIO
  #include <chemistry/qc/mbptr12/r12ia_mpiiofile.h>
#  endif
#endif

using namespace std;
using namespace sc;

inline int max(int a,int b) { return (a > b) ? a : b;}

/*-----------
  TwoBodyMOIntsTransform_ijxy
 -----------*/
static ClassDesc TwoBodyMOIntsTransform_ijxy_cd(
  typeid(TwoBodyMOIntsTransform_ijxy),"TwoBodyMOIntsTransform_ijxy",1,"public TwoBodyMOIntsTransform",
  0, 0, create<TwoBodyMOIntsTransform_ijxy>);

TwoBodyMOIntsTransform_ijxy::TwoBodyMOIntsTransform_ijxy(const std::string& name, const Ref<MOIntsTransformFactory>& factory,
                                                         const Ref<TwoBodyIntDescr>& tbintdescr,
                                                         const Ref<MOIndexSpace>& space1, const Ref<MOIndexSpace>& space2,
                                                         const Ref<MOIndexSpace>& space3, const Ref<MOIndexSpace>& space4) :
  TwoBodyMOIntsTransform(name,factory,tbintdescr,space1,space2,space3,space4)
{
  init_vars();
}

TwoBodyMOIntsTransform_ijxy::TwoBodyMOIntsTransform_ijxy(StateIn& si) : TwoBodyMOIntsTransform(si)
{
  init_vars();
}

TwoBodyMOIntsTransform_ijxy::~TwoBodyMOIntsTransform_ijxy()
{
}

void
TwoBodyMOIntsTransform_ijxy::save_data_state(StateOut& so)
{
  TwoBodyMOIntsTransform::save_data_state(so);
}

//////////////////////////////////////////////////////
// Compute required (dynamic) memory
// for a given batch size of the transformation
//
// Only arrays allocated before exiting the loop over
// i-batches are included here  - only these arrays
// affect the batch size.
//////////////////////////////////////////////////////
distsize_t
TwoBodyMOIntsTransform_ijxy::compute_transform_dynamic_memory_(int ni) const
{
  int nproc = msg_->n();
  int nthread = thr_->nthread();

  int rank2 = space2_->rank();
  int nbasis2 = space2_->basis()->nbasis();
  int nfuncmax3 = space3_->basis()->max_nfunction_in_shell();
  int nfuncmax4 = space4_->basis()->max_nfunction_in_shell();
  int rank3 = space3_->rank();
  int nbasis4 = space4_->basis()->nbasis();

  // compute nij as nij on node 0, since nij on node 0 is >= nij on other nodes
  int nij = compute_nij(ni, rank2, nproc, 0);

  distsize_t memsize = sizeof(double)*(num_te_types()*((distsize_t)nthread * ni * nbasis2 * nfuncmax3 * nfuncmax4 // iqrs
						     + (distsize_t)ni * rank2 * nfuncmax3 * nfuncmax4  // ijrs
						     + (distsize_t)nij * rank3 * nbasis4 // ijxs - buffer of 3 q.t. and higher
						     // transformed integrals
						     )
				       + (distsize_t)rank3 * nbasis4 // xs or xy
				       );

  return memsize;
}

const size_t
TwoBodyMOIntsTransform_ijxy::memgrp_blksize() const
{
  const int nbasis3 = space3_->basis()->nbasis();
  const int rank3 = space3_->rank();
  const int dim3 = (nbasis3 > rank3) ? nbasis3 : rank3;
  const int nbasis4 = space4_->basis()->nbasis();
  const int rank4 = space4_->rank();
  const int dim4 = (nbasis4 > rank4) ? nbasis4 : rank4;
  return dim3*dim4*sizeof(double);
}

void
TwoBodyMOIntsTransform_ijxy::init_acc()
{
  if (ints_acc_.nonnull())
    return;

  int nij = compute_nij(batchsize_, space2_->rank(), msg_->n(), msg_->me());

  alloc_mem((size_t)num_te_types()*nij*memgrp_blksize());

  switch (ints_method_) {

  case MOIntsTransformFactory::StoreMethod::mem_only:
    if (npass_ > 1)
      throw std::runtime_error("TwoBodyMOIntsTransform_ijxy::init_acc() -- cannot use MemoryGrp-based accumulator in multi-pass transformations");
#if HAVE_R12IA_MEMGRP
    ints_acc_ = new R12IntsAcc_MemoryGrp(mem_, num_te_types(), space1_->rank(), space2_->rank(), space3_->rank(), space4_->rank());  // Hack to avoid using nfzc and nocc
#else
    assert(false);
#endif
    break;

  case MOIntsTransformFactory::StoreMethod::mem_posix:
    if (npass_ == 1) {
#if HAVE_R12IA_MEMGRP
      ints_acc_ = new R12IntsAcc_MemoryGrp(mem_, num_te_types(), space1_->rank(), space2_->rank(), space3_->rank(), space4_->rank());
#else
      assert(false);
#endif
      break;
    }
    // else use the next case
      
  case MOIntsTransformFactory::StoreMethod::posix:
    ints_acc_ = new R12IntsAcc_Node0File((file_prefix_+"."+name_).c_str(), num_te_types(),
                                         space1_->rank(), space2_->rank(), space3_->rank(), space4_->rank());
    break;

#if HAVE_MPIIO
  case MOIntsTransformFactory::StoreMethod::mem_mpi:
    if (npass_ == 1) {
#if HAVE_R12IA_MEMGRP
      ints_acc_ = new R12IntsAcc_MemoryGrp(mem_, num_te_types(), space1_->rank(), space2_->rank(), space3_->rank(), space4_->rank());
#else
      assert(false);
#endif
      break;
    }
    // else use the next case

  case MOIntsTransformFactory::StoreMethod::mpi:
#if HAVE_R12IA_MPIIO
    ints_acc_ = new R12IntsAcc_MPIIOFile_Ind(mem_, (file_prefix_+"."+name_).c_str(), num_te_types(),
                                             space1_->rank(), space2_->rank(), space3_->rank(), space4_->rank());
#else
    assert(false);
#endif
    break;
#endif
  
  default:
    throw std::runtime_error("TwoBodyMOIntsTransform_ijxy::init_acc() -- invalid integrals store method");
  }
}

void
TwoBodyMOIntsTransform_ijxy::check_int_symm(double threshold) throw (ProgrammingError)
{
  Ref<R12IntsAcc> iacc = ints_acc();
  iacc->activate();

  int num_te_types = iacc->num_te_types();
  int ni = iacc->ni();
  int nj = iacc->nj();
  int nx = iacc->nx();
  int ny = iacc->ny();
  vector<unsigned int> isyms = space1_->mosym();
  vector<unsigned int> jsyms = space2_->mosym();
  vector<unsigned int> xsyms = space3_->mosym();
  vector<unsigned int> ysyms = space4_->mosym();
  
  int me = msg_->me();
  vector<int> twi_map;
  int ntasks_with_ints = iacc->tasks_with_access(twi_map);
  if (!iacc->has_access(me))
    return;
  
  int ij=0;
  for(int i=0; i<ni; i++) {
    int isym = isyms[i];
    for(int j=0; j<nj; j++, ij++) {
      int jsym = jsyms[j];
      if (ij%ntasks_with_ints != twi_map[me])
        continue;

      for(int t=0; t<num_te_types; t++) {
        const double* ints = iacc->retrieve_pair_block(i,j,static_cast<R12IntsAcc::tbint_type>(t));
        int xy=0;
        for(int x=0; x<nx; x++) {
          int xsym = xsyms[x];
          for(int y=0; y<ny; y++, xy++) {
            int ysym = ysyms[y];
            if ( (isym^jsym^xsym^ysym) != 0 && fabs(ints[xy]) > threshold) {
              ExEnv::outn() << scprintf("Integral type=%d i=%d j=%d x=%d y=%d should be zero\n",t,i,j,x,y);
              throw ProgrammingError("TwoBodyMOIntsTransform_ijxy::check_int_symm() -- nonzero nonsymmetric integrals are detected",
                                     __FILE__, __LINE__);
            }
          }
        }
        iacc->release_pair_block(i,j,static_cast<R12IntsAcc::tbint_type>(t));
      }
    }
  }
}

/*void
TwoBodyMOIntsTransform_ijxy::compute()
{
  init_acc();
}*/

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
