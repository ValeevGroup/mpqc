//
// r12ia_memgrp.cc
//
// Copyright (C) 2002 Edward Valeev
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <stdexcept>
#include <stdlib.h>
#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <chemistry/qc/mbptr12/r12ia_memgrp.h>

using namespace std;
using namespace sc;

///////////////////////////////////////////////////////////////

static ClassDesc R12IntsAcc_MemoryGrp_cd(
  typeid(R12IntsAcc_MemoryGrp),"R12IntsAcc_MemoryGrp",1,"public R12IntsAcc",
  0, 0, create<R12IntsAcc_MemoryGrp>);

R12IntsAcc_MemoryGrp::R12IntsAcc_MemoryGrp(Ref<MemoryGrp>& mem, int num_te_types, int nbasis1, int nbasis2, int nocc, int nfzc) :
  R12IntsAcc(num_te_types, nbasis1, nbasis2, nocc, nfzc)
{
  mem_ = mem;
  
  init();
}

R12IntsAcc_MemoryGrp::R12IntsAcc_MemoryGrp(StateIn& si) : R12IntsAcc(si)
{
  mem_ = MemoryGrp::get_default_memorygrp();
  
  init();
}

R12IntsAcc_MemoryGrp::~R12IntsAcc_MemoryGrp()
{
  for(int i=0;i<nocc_act_;i++)
    for(int j=0;j<nocc_act_;j++)
      if (!is_local(i,j)) {
	int ij = ij_index(i,j);
	for(int oper_type=0; oper_type<num_te_types(); oper_type++)
	  if (pairblk_[ij].ints_[oper_type] != NULL) {
	    ExEnv::outn() << indent << mem_->me() << ": i = " << i << " j = " << j << " oper_type = " << oper_type << endl;
	    throw std::runtime_error("Logic error: R12IntsAcc_MemoryGrp::~ : some nonlocal blocks have not been released!");
	  }
      }
  delete[] pairblk_;
}

void
R12IntsAcc_MemoryGrp::save_data_state(StateOut& so)
{
  R12IntsAcc::save_data_state(so);
}

void
R12IntsAcc_MemoryGrp::init()
{
  nproc_ = mem_->n();

  // Now do some extra work to figure layout of data in MemoryGrp
  // Compute global offsets to each processor's data
  int i,j,ij;
  pairblk_ = new struct PairBlkInfo[nocc_act_*nocc_act_];
  for(i=0,ij=0;i<nocc_act_;i++)
    for(j=0;j<nocc_act_;j++,ij++) {
      pairblk_[ij].ints_[eri] = NULL;
      pairblk_[ij].ints_[r12] = NULL;
      pairblk_[ij].ints_[r12t1] = NULL;
      pairblk_[ij].ints_[r12t2] = NULL;
      pairblk_[ij].refcount_[eri] = 0;
      pairblk_[ij].refcount_[r12] = 0;
      pairblk_[ij].refcount_[r12t1] = 0;
      pairblk_[ij].refcount_[r12t2] = 0;
      int local_ij_index = ij_index(i,j)/nproc_;
      pairblk_[ij].offset_ = (distsize_t)local_ij_index*blocksize_ + mem_->offset(ij_proc(i,j));
    }
}

void
R12IntsAcc_MemoryGrp::store_memorygrp(Ref<MemoryGrp>& mem, int ni)
{
  if (committed_) {
    ExEnv::out0() << "R12IntsAcc_MemoryGrp::store_memorygrp(mem,ni) called after all data has been committed" << endl;
    abort();
  }
  // mem must be the same as mem_
  else if (mem_ != mem) {
    ExEnv::out0() << "R12IntsAcc_MemoryGrp::store_memorygrp(mem,ni) called with invalid argument:" << endl <<
      "mem != R12IntsAcc_MemoryGrp::mem_" << endl;
    abort();
  }
  else if (ni != nocc_act_) {
    ExEnv::out0() << "R12IntsAcc_MemoryGrp::store_memorygrp(mem,ni) called with invalid argument:" << endl <<
      "ni != R12IntsAcc_MemoryGrp::nocc_act_" << endl;
    abort();
  }
  else
    for (int i=0; i<nocc_act_; i++)
      for (int j=0; j<nocc_act_; j++)
	if (is_local(i,j)) {
	  int local_ij_index = ij_index(i,j)/nproc_;
	  double *integral_ij_offset = (double *)mem_->localdata() + nbasis__2_*num_te_types()*local_ij_index;
	  store_pair_block(i,j,integral_ij_offset);
	}
}

void
R12IntsAcc_MemoryGrp::store_pair_block(int i, int j, double *ints)
{
  // For now store blocks local to this node ONLY
  if (is_local(i,j)) {
    int ij = ij_index(i,j);
    pairblk_[ij].ints_[eri] = ints;
    pairblk_[ij].ints_[r12] = ints + nbasis__2_;
    pairblk_[ij].ints_[r12t1] = ints + 2*nbasis__2_;
    pairblk_[ij].ints_[r12t2] = ints + 3*nbasis__2_;
  }
}

void
R12IntsAcc_MemoryGrp::deactivate()
{
  mem_->sync();
  mem_->set_localsize(0);
}
    
double *
R12IntsAcc_MemoryGrp::retrieve_pair_block(int i, int j, tbint_type oper_type)
{
  int ij = ij_index(i,j);
  struct PairBlkInfo *pb = &pairblk_[ij];
  if (!is_local(i,j) && pb->ints_[oper_type] == 0) {
    pb->ints_[oper_type] = (double *) mem_->obtain_readonly(pb->offset_ + (distsize_t)oper_type*blksize_, blksize_);
  }
  pb->refcount_[oper_type] += 1;
  return pb->ints_[oper_type];
}

void
R12IntsAcc_MemoryGrp::release_pair_block(int i, int j, tbint_type oper_type)
{
  int ij = ij_index(i,j);
  struct PairBlkInfo *pb = &pairblk_[ij];
  if (pb->refcount_[oper_type] <= 0) {
    ExEnv::outn() << indent << mem_->me() << ":refcount=0: i = " << i << " j = " << j << " tbint_type = " << oper_type << endl;
    throw std::runtime_error("Logic error: R12IntsAcc_MemoryGrp::release_pair_block: refcount is already zero!");
  }
  if (!is_local(i,j) && pb->ints_[oper_type] != NULL && pb->refcount_[oper_type] == 1) {
    mem_->release_readonly(pb->ints_[oper_type],pb->offset_+ oper_type*blksize_,blksize_);
    pb->ints_[oper_type] = NULL;
  }
  pb->refcount_[oper_type] -= 1;
}

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
