//
// r12ia_mpiiofile.cc
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
#include <chemistry/qc/mbptr12/r12ia_mpiiofile.h>

using namespace std;
using namespace sc;

///////////////////////////////////////////////////////////////

R12IntsAcc_MPIIOFile::R12IntsAcc_MPIIOFile(Ref<MemoryGrp>& mem, const char* filename, int nte_types, int nbasis1, int nbasis2,
					   int nocc, int nfzc, bool restart) :
    R12IntsAcc(nte_types, nbasis1, nbasis2, nocc, nfzc), datafile_(MPI_FILE_NULL)
{
  mem_ = mem;
  int errcod;
  errcod = MPI_Comm_size(MPI_COMM_WORLD, &nproc_);
  nints_per_block_ = nbasis__2_*num_te_types();

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
      pairblk_[ij].offset_ = (MPI_Offset)ij*blocksize_;
    }

  // Create the file
  icounter_ = 0;
  filename_ = strdup(filename);

  //
  // Try opening the file
  //
  int amode;
  if (!restart)
    amode = MPI_MODE_CREATE | MPI_MODE_DELETE_ON_CLOSE | MPI_MODE_WRONLY;
  else
    amode = MPI_MODE_RDONLY;

  errcod = MPI_File_open(MPI_COMM_WORLD, filename_, amode,
                         MPI_INFO_NULL, &datafile_);
  check_error_code_(errcod);
  errcod = MPI_File_close(&datafile_);
  check_error_code_(errcod);
}

R12IntsAcc_MPIIOFile::~R12IntsAcc_MPIIOFile()
{
  for(int i=0;i<nocc_act_;i++)
    for(int j=0;j<nocc_act_;j++) {
      int ij = ij_index(i,j);
      for(int oper_type=0; oper_type<num_te_types(); oper_type++)
	if (pairblk_[ij].ints_[oper_type] != NULL) {
	  ExEnv::outn() << indent << mem_->me() << ": i = " << i << " j = " << j << " oper_type = " << oper_type << endl;
	  throw std::runtime_error("Logic error: R12IntsAcc_MPIIOFile::~ : some nonlocal blocks have not been released!");
	}
    }
  delete[] pairblk_;
  free(filename_);
}


void
R12IntsAcc_MPIIOFile::check_error_code_(int errcod) const
{

  if (errcod != MPI_SUCCESS) {

    int errclass;
    MPI_Error_class(errcod, &errclass);

    switch (errclass) {
    default:
      char* errstr = new char[MPI_MAX_ERROR_STRING];
      int errstrlen;
      MPI_Error_string(errcod, errstr, &errstrlen);
      ExEnv::out0() << "R12IntsAcc_MPIIOFile::R12IntsAcc_MPIIOFile -- MPI-I/O error: " << errstr 
                    << " on file " << filename_ << endl;
      delete[] errstr;
      throw std::runtime_error("R12IntsAcc_MPIIOFile::R12IntsAcc_MPIIOFile -- MPI-I/O error");
    }
  }

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
  int errcod = MPI_File_open(MPI_COMM_WORLD, filename_, MPI_MODE_RDONLY | MPI_MODE_DELETE_ON_CLOSE, MPI_INFO_NULL, &datafile_);
  check_error_code_(errcod);
}

void
R12IntsAcc_MPIIOFile::deactivate()
{
  mem_->activate();
  int errcod = MPI_File_close(&datafile_);
  check_error_code_(errcod);
}

void
R12IntsAcc_MPIIOFile::release_pair_block(int i, int j, tbint_type oper_type)
{
  int ij = ij_index(i,j);
  struct PairBlkInfo *pb = &pairblk_[ij];
  if (pb->refcount_[oper_type] <= 0) {
    ExEnv::outn() << indent << mem_->me() << ":refcount=0: i = " << i << " j = " << j << " tbint_type = " << oper_type << endl;
    throw std::runtime_error("Logic error: R12IntsAcc_MPIIOFile::release_pair_block: refcount is already zero!");
  }
  if (pb->ints_[oper_type] != NULL && pb->refcount_[oper_type] == 1) {
    delete[] pb->ints_[oper_type];
    pb->ints_[oper_type] = NULL;
  }
  pb->refcount_[oper_type] -= 1;
}

///////////////////////////////////////////////////////////////

R12IntsAcc_MPIIOFile_Ind::R12IntsAcc_MPIIOFile_Ind(Ref<MemoryGrp>& mem, const char* filename, int num_te_types, int nbasis1, int nbasis2,
						   int nocc, int nfzc, bool restart) :
  R12IntsAcc_MPIIOFile(mem,filename,num_te_types,nbasis1,nbasis2,nocc,nfzc,restart)
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
    check_error_code_(errcod);
    
    for (int i=0; i<ni; i++)
      for (int j=0; j<nocc_act_; j++) {
	int ij = ij_index(i,j);
	int proc = ij%nproc;
        if (proc != me) continue;
	int local_ij_index = ij/nproc;
        double *data = (double *) mem->localdata() + nbasis__2_*num_te_types()*local_ij_index;
        
        int IJ = ij_index(i+icounter_,j);
        int errcod = MPI_File_seek(datafile_, pairblk_[IJ].offset_, MPI_SEEK_SET);
        check_error_code_(errcod);
        MPI_Status status;
        errcod = MPI_File_write(datafile_, (void *)data, nints_per_block_, MPI_DOUBLE, &status);
        check_error_code_(errcod);
      }
    // Close the file and update the i counter
    errcod = MPI_File_close(&datafile_);
    check_error_code_(errcod);
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
    int errcod = MPI_File_seek(datafile_, offset, MPI_SEEK_SET);
    check_error_code_(errcod);
    double *buffer = new double[nbasis__2_];
    MPI_Status status;
    errcod = MPI_File_read(datafile_, (void *)buffer, nbasis__2_, MPI_DOUBLE, &status);
    check_error_code_(errcod);
    pb->ints_[oper_type] = buffer;
  }
  pb->refcount_[oper_type] += 1;
  return pb->ints_[oper_type];
}

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
