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
#include <util/misc/scexception.h>
#include <chemistry/qc/mbptr12/r12ia_memgrp.h>

using namespace std;
using namespace sc;

///////////////////////////////////////////////////////////////

static ClassDesc R12IntsAcc_MemoryGrp_cd(
  typeid(R12IntsAcc_MemoryGrp),"R12IntsAcc_MemoryGrp",1,"public R12IntsAcc",
  0, 0, create<R12IntsAcc_MemoryGrp>);

R12IntsAcc_MemoryGrp::R12IntsAcc_MemoryGrp(Ref<MemoryGrp>& mem, int num_te_types, int ni, int nj, int nx, int ny) :
  R12IntsAcc(num_te_types, ni, nj, nx, ny), blksize_memgrp_(blksize_)
{
  mem_ = mem;
  
  init();
}

R12IntsAcc_MemoryGrp::R12IntsAcc_MemoryGrp(StateIn& si) : R12IntsAcc(si)
{
  mem_ = MemoryGrp::get_default_memorygrp();
  blksize_memgrp_ = blksize_;
  
  init();
}

R12IntsAcc_MemoryGrp::~R12IntsAcc_MemoryGrp()
{
  for(int i=0;i<ni_;i++)
    for(int j=0;j<nj_;j++)
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
  pairblk_ = new struct PairBlkInfo[ni_*nj_];
  for(i=0,ij=0;i<ni_;i++)
    for(j=0;j<nj_;j++,ij++) {
      for(int type=0; type<num_te_types(); type++) {
        pairblk_[ij].ints_[type] = NULL;
        pairblk_[ij].refcount_[type] = 0;
        }
      int local_ij_index = ij_index(i,j)/nproc_;
      pairblk_[ij].offset_ = (distsize_t)local_ij_index*blksize_memgrp_*num_te_types() +
                             mem_->offset(ij_proc(i,j));
      }
}

void
R12IntsAcc_MemoryGrp::store_memorygrp(Ref<MemoryGrp>& mem, int ni, const size_t blksize)
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
  else if (ni != ni_) {
    ExEnv::out0() << "R12IntsAcc_MemoryGrp::store_memorygrp(mem,ni) called with invalid argument:" << endl <<
      "ni != R12IntsAcc_MemoryGrp::ni_" << endl;
    abort();
  }
  else {
    if (blksize != 0 && blksize != blksize_memgrp_) {
      blksize_memgrp_ = blksize;
      init();
    }
    
    for (int i=0; i<ni_; i++)
      for (int j=0; j<nj_; j++)
	if (is_local(i,j)) {
	  int local_ij_index = ij_index(i,j)/nproc_;
	  double *integral_ij_offset = (double *) ((size_t)mem_->localdata() +
                                            blksize_memgrp_*num_te_types()*local_ij_index);
	  store_pair_block(i,j,integral_ij_offset);
	}
  }

  inc_next_orbital(ni);
}

void
R12IntsAcc_MemoryGrp::store_pair_block(int i, int j, double *ints)
{
  // For now store blocks local to this node ONLY
  if (is_local(i,j)) {
    int ij = ij_index(i,j);
    pairblk_[ij].ints_[0] = ints;
    const size_t blksize = blksize_memgrp_/sizeof(double);
    for(int type=1; type<num_te_types(); type++)
      pairblk_[ij].ints_[type] = pairblk_[ij].ints_[type-1] + blksize;
  }
}

void
R12IntsAcc_MemoryGrp::deactivate()
{
  R12IntsAcc::deactivate();
  mem_->sync();
  mem_->set_localsize(0);
  // so that this accumator cannot be activated again
  committed_ = false;
}
    
double *
R12IntsAcc_MemoryGrp::retrieve_pair_block(int i, int j, tbint_type oper_type)
{
  int ij = ij_index(i,j);
  struct PairBlkInfo *pb = &pairblk_[ij];
  if (!is_local(i,j) && pb->ints_[oper_type] == 0) {
    pb->ints_[oper_type] = (double *) mem_->obtain_readonly(pb->offset_ + (distsize_t)oper_type*blksize_memgrp_, blksize_);
  }
  pb->refcount_[oper_type] += 1;
  if (classdebug() > 0)
    ExEnv::outn() << indent << mem_->me() << ":refcount=" << pb->refcount_[oper_type]
                  << ": i = " << i << " j = " << j << " tbint_type = " << oper_type << endl;
  return pb->ints_[oper_type];
}

void
R12IntsAcc_MemoryGrp::release_pair_block(int i, int j, tbint_type oper_type)
{
  int ij = ij_index(i,j);
  struct PairBlkInfo *pb = &pairblk_[ij];
  if (pb->refcount_[oper_type] <= 0) {
    ExEnv::outn() << indent << mem_->me() << ":refcount=0: i = " << i << " j = " << j << " tbint_type = " << oper_type << endl;
    throw ProgrammingError("Logic error: R12IntsAcc_MemoryGrp::release_pair_block: refcount is already zero!",__FILE__,__LINE__);
  }
  if (!is_local(i,j) && pb->ints_[oper_type] != NULL && pb->refcount_[oper_type] == 1) {
    mem_->release_readonly(pb->ints_[oper_type],pb->offset_+ oper_type*blksize_memgrp_,blksize_);
    pb->ints_[oper_type] = NULL;
  }
  pb->refcount_[oper_type] -= 1;
  if (classdebug() > 0)
    ExEnv::outn() << indent << mem_->me() << ":refcount=" << pb->refcount_[oper_type]
                  << ": i = " << i << " j = " << j << " tbint_type = " << oper_type << endl;
}

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
