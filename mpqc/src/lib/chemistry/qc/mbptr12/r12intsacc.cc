//
// r12intsacc.cc
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

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <chemistry/qc/mbptr12/r12intsacc.h>

using namespace std;
using namespace sc;

/*--------------------------------
  R12IntsAcc
 --------------------------------*/

R12IntsAcc::R12IntsAcc(int nbasis, int noso, int nocc, int nfzc)
{
  noso_ = noso;
  nocc_ = nocc;
  nfzc_ = nfzc;
  nocc_act_ = nocc_ - nfzc_;
  nbasis_ = nbasis;
  nbasis__2_ = nbasis_*nbasis_;
  blksize_ = nbasis__2_*sizeof(double);
  blocksize_ = blksize_*num_te_types_;
  committed_ = false;
}

R12IntsAcc::~R12IntsAcc()
{
}


///////////////////////////////////////////////////////////////

R12IntsAcc_MemoryGrp::R12IntsAcc_MemoryGrp(Ref<MemoryGrp>& mem, int nbasis, int noso, int nocc, int nfzc) :
  R12IntsAcc(nbasis, noso, nocc, nfzc)
{
  mem_ = mem;
  nproc_ = mem->n();
  
  // Now do some extra work to figure layout of data in MemoryGrp
  // Compute global offsets to each processor's data
  int i,j,ij;
  pairblk_ = new struct PairBlkInfo[nocc_act_*nocc_act_];
  for(i=0,ij=0;i<nocc_act_;i++)
    for(j=0;j<nocc_act_;j++,ij++) {
      pairblk_[ij].ints_[eri] = NULL;
      pairblk_[ij].ints_[r12] = NULL;
      pairblk_[ij].ints_[r12t1] = NULL;
      int local_ij_index = ij_index(i,j)/nproc_;
      pairblk_[ij].offset_ = (distsize_t)local_ij_index*blocksize_ + mem_->offset(ij_proc(i,j));
    }
}

R12IntsAcc_MemoryGrp::~R12IntsAcc_MemoryGrp()
{
  for(int i=0;i<nocc_act_;i++)
    for(int j=0;j<nocc_act_;j++) {
      release_pair_block(i,j,eri);
      release_pair_block(i,j,r12);
      release_pair_block(i,j,r12t1);
    }
  delete[] pairblk_;
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
	  double *integral_ij_offset = (double *)mem_->localdata() + nbasis__2_*num_te_types_*local_ij_index;
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
  // Can retrieve both local ...
  if (is_local(i,j)) {
    switch(oper_type) {
      // Locally subblocks are store next to each other
    case eri:
      return pb->ints_[eri];
    case r12:
      return pb->ints_[r12];
    case r12t1:
      return pb->ints_[r12t1];
    default:
      return 0;
    }
  }
  // ... and remote blocks
  else {
    if (pb->ints_[oper_type] == 0)
      pb->ints_[oper_type] = (double *) mem_->obtain_readonly(pb->offset_ + (distsize_t)oper_type*blksize_, blksize_);
    return pb->ints_[oper_type];
  }
}

void
R12IntsAcc_MemoryGrp::release_pair_block(int i, int j, tbint_type oper_type)
{
  if (is_local(i,j)) {
    // do nothing
  }
  else {
    int ij = ij_index(i,j);
    struct PairBlkInfo *pb = &pairblk_[ij];
    if (pb->ints_[oper_type] != NULL) {
      mem_->release_readonly(pb->ints_[oper_type],pb->offset_+ oper_type*blksize_,blksize_);
      pb->ints_[oper_type] = NULL;
    }
  }
}


///////////////////////////////////////////////////////////////

R12IntsAcc_Node0File::R12IntsAcc_Node0File(Ref<MemoryGrp>& mem, char* filename, int nbasis, int noso, int nocc, int nfzc, bool restart) :
  R12IntsAcc(nbasis, noso, nocc, nfzc)
{
  mem_ = mem;

  pairblk_ = new struct PairBlkInfo[nocc_act_*nocc_act_];
  int i, j, ij;
  for(i=0,ij=0;i<nocc_act_;i++)
    for(j=0;j<nocc_act_;j++,ij++) {
      pairblk_[ij].ints_[eri] = NULL;
      pairblk_[ij].ints_[r12] = NULL;
      pairblk_[ij].ints_[r12t1] = NULL;
      pairblk_[ij].offset_ = (off_t)ij*blocksize_;
    }

  // Create the file
  icounter_ = 0;
  filename_ = strdup(filename);
  if (!restart) {
    datafile_ = creat(filename_,0644);
    close(datafile_);
  }
}

R12IntsAcc_Node0File::~R12IntsAcc_Node0File()
{
  for(int i=0;i<nocc_act_;i++)
    for(int j=0;j<nocc_act_;j++) {
      release_pair_block(i,j,eri);
      release_pair_block(i,j,r12);
      release_pair_block(i,j,r12t1);
    }
  delete[] pairblk_;
  delete[] filename_;
      
  // Destroy the file
  unlink(filename_);
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
  else if (icounter_ + ni > nocc_act_) {
    ExEnv::out0() << "R12IntsAcc_Node0File::store_memorygrp(mem,ni) called with invalid argument:" << endl <<
      "ni+icounter_ > R12IntsAcc_Node0File::nocc_act_" << endl;
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
	  data = (double *) mem->localdata() + nbasis__2_*num_te_types_*local_ij_index;
	  write(datafile_, data, blocksize_);
	}
      }
    // Close the file and update the i counter
    close(datafile_);
    icounter_ += ni;
  }
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
    return pb->ints_[oper_type];
  }
  else {
    ExEnv::out0() << "R12IntsAcc_Node0File::retrieve_pair_block() called on node other than 0" << endl;
    abort();
  }
  return 0;
}

void
R12IntsAcc_Node0File::release_pair_block(int i, int j, tbint_type oper_type)
{
  if (is_avail(i,j)) {
    int ij = ij_index(i,j);
    struct PairBlkInfo *pb = &pairblk_[ij];
    if (pb->ints_[oper_type] != NULL) {
      delete[] pb->ints_[oper_type];
      pb->ints_[oper_type] = NULL;
    }
  }
}

///////////////////////////////////////////////////////////////

R12IntsAcc_MPIIOFile::R12IntsAcc_MPIIOFile(Ref<MemoryGrp>& mem, char* filename, int nbasis, int noso, int nocc, int nfzc, bool restart) :
    R12IntsAcc(nbasis, noso, nocc, nfzc), datafile_(MPI_FILE_NULL)
{
  mem_ = mem;
  int errcod;
  errcod = MPI_Comm_size(MPI_COMM_WORLD, &nproc_);
  nints_per_block_ = nbasis__2_*num_te_types_;

  pairblk_ = new struct PairBlkInfo[nocc_act_*nocc_act_];
  int i, j, ij;
  for(i=0,ij=0;i<nocc_act_;i++)
    for(j=0;j<nocc_act_;j++,ij++) {
      pairblk_[ij].ints_[eri] = NULL;
      pairblk_[ij].ints_[r12] = NULL;
      pairblk_[ij].ints_[r12t1] = NULL;
      pairblk_[ij].offset_ = (MPI_Offset)ij*blocksize_;
    }

  // Create the file
  icounter_ = 0;
  filename_ = strdup(filename);
  if (!restart) {
    MPI_File_open(MPI_COMM_WORLD, filename_, MPI_MODE_CREATE | MPI_MODE_DELETE_ON_CLOSE, MPI_INFO_NULL, &datafile_);
    MPI_File_close(&datafile_);
  }
}

R12IntsAcc_MPIIOFile::~R12IntsAcc_MPIIOFile()
{
  for(int i=0;i<nocc_act_;i++)
    for(int j=0;j<nocc_act_;j++) {
      release_pair_block(i,j,eri);
      release_pair_block(i,j,r12);
      release_pair_block(i,j,r12t1);
    }
  delete[] pairblk_;
  delete[] filename_;
}

void
R12IntsAcc_MPIIOFile::store_pair_block(int i, int j, double *ints)
{
  ExEnv::err0() << "R12IntsAcc_MPIIOFile::store_pair_block() called: error" << endl;
  abort();
}

void
R12IntsAcc_MPIIOFile::commit()
{
  mem_->set_localsize(0);
  mem_->sync();
  mem_->deactivate();
  R12IntsAcc::commit();
  MPI_File_open(MPI_COMM_WORLD, filename_, MPI_MODE_RDONLY | MPI_MODE_DELETE_ON_CLOSE, MPI_INFO_NULL, &datafile_);
}

void
R12IntsAcc_MPIIOFile::deactivate()
{
  mem_->activate();
  MPI_File_close(&datafile_);
}

void
R12IntsAcc_MPIIOFile::release_pair_block(int i, int j, tbint_type oper_type)
{
  int ij = ij_index(i,j);
  struct PairBlkInfo *pb = &pairblk_[ij];
  if (pb->ints_[oper_type] != NULL) {
    delete[] pb->ints_[oper_type];
    pb->ints_[oper_type] = NULL;
  }
}

///////////////////////////////////////////////////////////////

R12IntsAcc_MPIIOFile_Ind::R12IntsAcc_MPIIOFile_Ind(Ref<MemoryGrp>& mem, char* filename, int nbasis, int noso, int nocc, int nfzc, bool restart) :
  R12IntsAcc_MPIIOFile(mem,filename,nbasis,noso,nocc,nfzc,restart)
{
}

R12IntsAcc_MPIIOFile_Ind::~R12IntsAcc_MPIIOFile_Ind()
{
}

void
R12IntsAcc_MPIIOFile_Ind::store_memorygrp(Ref<MemoryGrp>& mem, int ni)
{
  if (committed_) {
    ExEnv::out0() << "R12IntsAcc_MPIIOFile_Ind::store_memorygrp(mem,ni) called after all data has been committed" << endl;
    abort();
  }
  // mem must be the same as mem_
  else if (mem_ != mem) {
    ExEnv::out0() << "R12IntsAcc_MemoryGrp::store_memorygrp(mem,ni) called with invalid argument:" << endl <<
      "mem != R12IntsAcc_MemoryGrp::mem_" << endl;
    abort();
  }
  else if (ni > nocc_act_) {
    ExEnv::out0() << "R12IntsAcc_MPIIOFile_Ind::store_memorygrp(mem,ni) called with invalid argument:" << endl <<
      "ni > R12IntsAcc_MPIIOFile_Ind::nocc_act_" << endl;
    abort();
  }
  else if (icounter_ + ni > nocc_act_) {
    ExEnv::out0() << "R12IntsAcc_MPIIOFile_Ind::store_memorygrp(mem,ni) called with invalid argument:" << endl <<
      "ni+icounter_ > R12IntsAcc_MPIIOFile_Ind::nocc_act_" << endl;
    abort();
  }
  else {
    // Now do some extra work to figure layout of data in MemoryGrp
    // Compute global offsets to each processor's data
    int i,j,ij;
    int me = mem->me();
    int nproc = mem->n();

    // Append the data to the file
    int errcod = MPI_File_open(MPI_COMM_WORLD, filename_, MPI_MODE_CREATE | MPI_MODE_APPEND | MPI_MODE_WRONLY, MPI_INFO_NULL, &datafile_);
    
    for (int i=0; i<ni; i++)
      for (int j=0; j<nocc_act_; j++) {
	int ij = ij_index(i,j);
	int proc = ij%nproc;
        if (proc != me) continue;
	int local_ij_index = ij/nproc;
        double *data = (double *) mem->localdata() + nbasis__2_*num_te_types_*local_ij_index;
        
        int IJ = ij_index(i+icounter_,j);
        MPI_File_seek(datafile_, pairblk_[IJ].offset_, MPI_SEEK_SET);
        MPI_Status status;
        MPI_File_write(datafile_, (void *)data, nints_per_block_, MPI_DOUBLE, &status);
      }
    // Close the file and update the i counter
    MPI_File_close(&datafile_);
    icounter_ += ni;
  }
}

double *
R12IntsAcc_MPIIOFile_Ind::retrieve_pair_block(int i, int j, tbint_type oper_type)
{
  int ij = ij_index(i,j);
  struct PairBlkInfo *pb = &pairblk_[ij];
  // Always first check if it's already in memory
  if (pb->ints_[oper_type] == 0) {
    MPI_Offset offset = pb->offset_ + (MPI_Offset)oper_type*blksize_;
    MPI_File_seek(datafile_, offset, MPI_SEEK_SET);
    double *buffer = new double[nbasis__2_];
    MPI_Status status;
    MPI_File_read(datafile_, (void *)buffer, nbasis__2_, MPI_DOUBLE, &status);
    pb->ints_[oper_type] = buffer;
  }
  return pb->ints_[oper_type];
}

///////////////////////////////////////////////////////////////

R12IntsAcc_MPIIOFile_Coll::R12IntsAcc_MPIIOFile_Coll(Ref<MemoryGrp>& mem, char* filename, int nbasis, int noso, int nocc, int nfzc, bool restart) :
  R12IntsAcc_MPIIOFile(mem,filename,nbasis,noso,nocc,nfzc,restart)
{
  ExEnv::out0() << "R12IntsAcc_MPIIOFile_Coll::R12IntsAcc_MPIIOFile_Coll has been called but this class should not be used" << endl;
  abort();
}

R12IntsAcc_MPIIOFile_Coll::~R12IntsAcc_MPIIOFile_Coll()
{
}

void
R12IntsAcc_MPIIOFile_Coll::store_memorygrp(Ref<MemoryGrp>& mem, int ni)
{
  if (committed_) {
    ExEnv::out0() << "R12IntsAcc_MPIIOFile_Coll::store_memorygrp(mem,ni) called after all data has been committed" << endl;
    abort();
  }
  // mem must be the same as mem_
  else if (mem_ != mem) {
    ExEnv::out0() << "R12IntsAcc_MemoryGrp::store_memorygrp(mem,ni) called with invalid argument:" << endl <<
      "mem != R12IntsAcc_MemoryGrp::mem_" << endl;
    abort();
  }
  else if (ni > nocc_act_) {
    ExEnv::out0() << "R12IntsAcc_MPIIOFile_Coll::store_memorygrp(mem,ni) called with invalid argument:" << endl <<
      "ni > R12IntsAcc_MPIIOFile_Coll::nocc_act_" << endl;
    abort();
  }
  else if (icounter_ + ni > nocc_act_) {
    ExEnv::out0() << "R12IntsAcc_MPIIOFile_Coll::store_memorygrp(mem,ni) called with invalid argument:" << endl <<
      "ni+icounter_ > R12IntsAcc_MPIIOFile_Coll::nocc_act_" << endl;
    abort();
  }
  else {
    // Now do some extra work to figure layout of data in MemoryGrp
    // Compute global offsets to each processor's data
    int i,j,ij;
    int errcod;
    int me = mem->me();
    int nproc = mem->n();

    /*-----------------------------
      Append the data to the file
     -----------------------------*/

    // figure out how many blocks this node has
    int nij = ni*nocc_act_;
    int nij_local = nij/nproc;
    int remainder = nij%nproc;
    nij_local += (me < remainder) ? 1 : 0;
    int stride = nproc_ * nints_per_block_;
    errcod = MPI_Type_vector(nij_local, nints_per_block_, stride, MPI_DOUBLE, &filetype_);
    MPI_Type_commit(&filetype_);

    errcod = MPI_File_open(MPI_COMM_WORLD, filename_, MPI_MODE_CREATE | MPI_MODE_APPEND | MPI_MODE_WRONLY, MPI_INFO_NULL, &datafile_);
    MPI_Offset offset = (icounter_*nocc_act_ + me) * (MPI_Offset)blocksize_;
    errcod = MPI_File_set_view(datafile_, offset, MPI_DOUBLE, filetype_, "native", MPI_INFO_NULL);

    // Write the data out
    size_t nints_local = nij_local*(size_t)nints_per_block_;
    errcod = MPI_File_write_all(datafile_, (void *)mem->localdata(), nints_local, MPI_DOUBLE, MPI_STATUS_IGNORE);
    
    // Close the file and update the i counter
    MPI_File_close(&datafile_);
    icounter_ += ni;
    MPI_Type_free(&filetype_);
  }
}

double *
R12IntsAcc_MPIIOFile_Coll::retrieve_pair_block(int i, int j, tbint_type oper_type)
{
  int ij = ij_index(i,j);
  struct PairBlkInfo *pb = &pairblk_[ij];
  // Always first check if it's already in memory
  if (pb->ints_[oper_type] == 0) {
    MPI_Offset offset = pb->offset_ + oper_type*blksize_;
    MPI_File_seek(datafile_, offset, MPI_SEEK_SET);
    double *buffer = new double[nbasis__2_];
    MPI_Status status;
    MPI_File_read(datafile_, (void *)buffer, blksize_, MPI_BYTE, &status);
    pb->ints_[oper_type] = buffer;
  }
  return pb->ints_[oper_type];
}

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
