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
#include <chemistry/qc/mbptr12/r12ia_memgrp.h>
#include <chemistry/qc/mbptr12/r12ia_node0file.h>
#ifdef HAVE_MPIIO
  #include <chemistry/qc/mbptr12/r12ia_mpiiofile.h>
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
                                                         const Ref<MOIndexSpace>& space1, const Ref<MOIndexSpace>& space2,
                                                         const Ref<MOIndexSpace>& space3, const Ref<MOIndexSpace>& space4) :
  TwoBodyMOIntsTransform(name,factory,space1,space2,space3,space4)
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

  distsize_t memsize = sizeof(double)*(num_te_types_*((distsize_t)nthread * ni * nbasis2 * nfuncmax3 * nfuncmax4 // iqrs
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

  alloc_mem((size_t)num_te_types_*nij*memgrp_blksize());

  // R12IntsAcc cannot work yet in cases when i and j are different spaces
  if (space1_ != space2_)
    throw std::runtime_error("TwoBodyMOIntsTransform_ijxy::init_acc() -- space1_ and space2_ must be the same");

  switch (ints_method_) {

  case MOIntsTransformFactory::mem_only:
    if (npass_ > 1)
      throw std::runtime_error("TwoBodyMOIntsTransform_ijxy::init_acc() -- cannot use MemoryGrp-based accumulator in multi-pass transformations");
    ints_acc_ = new R12IntsAcc_MemoryGrp(mem_, num_te_types_, space3_->rank(), space4_->rank(), space1_->rank());  // Hack to avoid using nfzc and nocc
    break;

  case MOIntsTransformFactory::mem_posix:
    if (npass_ == 1) {
      ints_acc_ = new R12IntsAcc_MemoryGrp(mem_, num_te_types_, space3_->rank(), space4_->rank(), space1_->rank());
      break;
    }
    // else use the next case
      
  case MOIntsTransformFactory::posix:
    ints_acc_ = new R12IntsAcc_Node0File(mem_, (file_prefix_+"."+name_).c_str(), num_te_types_,
                                         space3_->rank(), space4_->rank(), space1_->rank());
    break;

#if HAVE_MPIIO
  case MOIntsTransformFactory::mem_mpi:
    if (npass_ == 1) {
      ints_acc_ = new R12IntsAcc_MemoryGrp(mem_, num_te_types_, space3_->rank(), space4_->rank(), space1_->rank());
      break;
    }
    // else use the next case

  case MOIntsTransformFactory::mpi:
    ints_acc_ = new R12IntsAcc_MPIIOFile_Ind(mem_, (file_prefix_+"."+name_).c_str(), num_te_types_,
                                             space3_->rank(), space4_->rank(), space1_->rank());
    break;
#endif
  
  default:
    throw std::runtime_error("TwoBodyMOIntsTransform_ijxy::init_acc() -- invalid integrals store method");
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
