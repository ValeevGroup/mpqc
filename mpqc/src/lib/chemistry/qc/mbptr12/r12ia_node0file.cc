//
// r12ia_node0file.cc
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
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <stdlib.h>
#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <chemistry/qc/mbptr12/r12ia_node0file.h>

using namespace std;
using namespace sc;

///////////////////////////////////////////////////////////////

static ClassDesc R12IntsAcc_Node0File_cd(
  typeid(R12IntsAcc_Node0File),"R12IntsAcc_Node0File",1,"public R12IntsAcc",
  0, 0, create<R12IntsAcc_Node0File>);

R12IntsAcc_Node0File::R12IntsAcc_Node0File(Ref<MemoryGrp>& mem, const char* filename, int num_te_types,
					   int nbasis1, int nbasis2, int nocc_act, bool restart) :
  R12IntsAcc(num_te_types, nbasis1, nbasis2, nocc_act)
{
  mem_ = mem;
  filename_ = strdup(filename);
  
  init(restart);
}

R12IntsAcc_Node0File::R12IntsAcc_Node0File(StateIn& si) :
  R12IntsAcc(si)
{
  mem_ = MemoryGrp::get_default_memorygrp();
  si.getstring(filename_);
  
  init(true);
}

R12IntsAcc_Node0File::~R12IntsAcc_Node0File()
{
  for(int i=0;i<nocc_act_;i++)
    for(int j=0;j<nocc_act_;j++) {
      if (!is_avail(i,j)) {
	int ij = ij_index(i,j);
	for(int oper_type=0; oper_type<num_te_types(); oper_type++)
	  if (pairblk_[ij].ints_[oper_type] != NULL) {
	    ExEnv::outn() << indent << mem_->me() << ": i = " << i << " j = " << j << " oper_type = " << oper_type << endl;
	    throw std::runtime_error("Logic error: R12IntsAcc_Node0File::~ : some nonlocal blocks have not been released!");
	  }
      }
    }
  delete[] pairblk_;
      
  // Destroy the file
  unlink(filename_);
  free(filename_);
}

void
R12IntsAcc_Node0File::save_data_state(StateOut& so)
{
  R12IntsAcc::save_data_state(so);
  so.putstring(filename_);
}

void
R12IntsAcc_Node0File::init(bool restart)
{
  pairblk_ = new struct PairBlkInfo[nocc_act_*nocc_act_];
  int i, j, ij;
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
      pairblk_[ij].offset_ = (off_t)ij*blocksize_;
    }

  // See if can open/create the file
  if (restart)
    datafile_ = open(filename_,O_WRONLY|O_APPEND,0644);
  else
    datafile_ = creat(filename_,0644);
  // Check if the file was opened correctly
  check_filedescr_();
  // If everything is fine close it and proceed
  close(datafile_);
}

void
R12IntsAcc_Node0File::check_filedescr_()
{
  // If the file was not opened correctly - throw an exception
  if (datafile_ == -1) {
    switch (errno) {
      case EACCES:
      throw std::runtime_error("R12IntsAcc_Node0File::R12IntsAcc_Node0File -- access to the requested file is not allowed");
      case ENOSPC:
      throw std::runtime_error("R12IntsAcc_Node0File::R12IntsAcc_Node0File -- no space left in the filesystem");

      default:
      throw std::runtime_error("R12IntsAcc_Node0File::R12IntsAcc_Node0File -- failed to open POSIX file on node 0");
    }
  }
}

void
R12IntsAcc_Node0File::store_memorygrp(Ref<MemoryGrp>& mem, int ni)
{
  if (committed_) {
    ExEnv::out0() << "R12IntsAcc_Node0File::store_memorygrp(mem,ni) called after all data has been committed" << endl;
    abort();
  }
  // mem must be the same as mem_
  else if (mem_ != mem) {
    ExEnv::out0() << "R12IntsAcc_MemoryGrp::store_memorygrp(mem,ni) called with invalid argument:" << endl <<
      "mem != R12IntsAcc_MemoryGrp::mem_" << endl;
    abort();
  }
  // Will store integrals on node 0
  else if (mem->me() != 0)
    return;
  else if (ni > nocc_act_) {
    ExEnv::out0() << "R12IntsAcc_Node0File::store_memorygrp(mem,ni) called with invalid argument:" << endl <<
      "ni > R12IntsAcc_Node0File::nocc_act_" << endl;
    abort();
  }
  else if (next_orbital() + ni > nocc_act_) {
    ExEnv::out0() << "R12IntsAcc_Node0File::store_memorygrp(mem,ni) called with invalid argument:" << endl <<
      "ni+next_orbital() > R12IntsAcc_Node0File::nocc_act_" << endl;
    abort();
  }
  else {
    // Now do some extra work to figure layout of data in MemoryGrp
    // Compute global offsets to each processor's data
    int i,j,ij;
    int me = mem->me();
    int nproc = mem->n();

    // Append the data to the file
    datafile_ = open(filename_,O_WRONLY|O_APPEND,0644);
    for (int i=0; i<ni; i++)
      for (int j=0; j<nocc_act_; j++) {
	double *data;
	int ij = ij_index(i,j);
	int proc = ij%nproc;
	int local_ij_index = ij/nproc;
	if (proc != me) {
	  distsize_t moffset = (distsize_t)local_ij_index*blocksize_ + mem->offset(proc);
	  data = (double *) mem->obtain_readonly(moffset, blocksize_);
	  write(datafile_, data, blocksize_);
	  mem->release_readonly(data, moffset, blocksize_);
	}
	else {
	  data = (double *) mem->localdata() + nbasis__2_*num_te_types()*local_ij_index;
	  write(datafile_, data, blocksize_);
	}
      }
    // Close the file and update the i counter
    close(datafile_);
  }

  inc_next_orbital(ni);
}

void
R12IntsAcc_Node0File::store_pair_block(int i, int j, double *ints)
{
  ExEnv::err0() << "R12IntsAcc_Node0File::store_pair_block() called: error" << endl;
  abort();
}

void
R12IntsAcc_Node0File::commit()
{
  mem_->set_localsize(0);
  mem_->sync();
  mem_->deactivate();
  R12IntsAcc::commit();
  datafile_ = open(filename_, O_RDONLY);
}

void
R12IntsAcc_Node0File::deactivate()
{
  mem_->activate();
  close(datafile_);
}

double *
R12IntsAcc_Node0File::retrieve_pair_block(int i, int j, tbint_type oper_type)
{
  // Can retrieve blocks on node 0 only
  if (is_avail(i,j)) {
    int ij = ij_index(i,j);
    struct PairBlkInfo *pb = &pairblk_[ij];
    // Always first check if it's already in memory
    if (pb->ints_[oper_type] == 0) {
      off_t offset = pb->offset_ + (off_t)oper_type*blksize_;
      lseek(datafile_,offset,SEEK_SET);
      pb->ints_[oper_type] = new double[nbasis__2_];
      read(datafile_,pb->ints_[oper_type],blksize_);
    }
    pb->refcount_[oper_type] += 1;
    return pb->ints_[oper_type];
  }
  else
    throw std::runtime_error("R12IntsAcc_Node0File::retrieve_pair_block() called on node other than 0");

  return 0;
}

void
R12IntsAcc_Node0File::release_pair_block(int i, int j, tbint_type oper_type)
{
  if (is_avail(i,j)) {
    int ij = ij_index(i,j);
    struct PairBlkInfo *pb = &pairblk_[ij];
    if (pb->refcount_[oper_type] <= 0) {
      ExEnv::outn() << indent << mem_->me() << ":refcount=0: i = " << i << " j = " << j << " tbint_type = " << oper_type << endl;
      throw std::runtime_error("Logic error: R12IntsAcc_Node0File::release_pair_block: refcount is already zero!");
    }
    if (pb->ints_[oper_type] != NULL && pb->refcount_[oper_type] == 1) {
      delete[] pb->ints_[oper_type];
      pb->ints_[oper_type] = NULL;
    }
    pb->refcount_[oper_type] -= 1;
  }
}


// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
