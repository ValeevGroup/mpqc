//
// r12ia_mpiiofile.cc
//
// Copyright (C) 2002 Edward Valeev
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <cassert>
#include <stdexcept>
#include <stdlib.h>
#include <util/misc/string.h>
#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <util/class/scexception.h>
#include <chemistry/qc/mbptr12/r12ia_mpiiofile.h>

using namespace std;
using namespace sc;

///////////////////////////////////////////////////////////////

static ClassDesc R12IntsAcc_MPIIOFile_cd(
  typeid(R12IntsAcc_MPIIOFile),"R12IntsAcc_MPIIOFile",1,"public R12IntsAcc",
  0, 0, 0);

R12IntsAcc_MPIIOFile::R12IntsAcc_MPIIOFile(const char* filename, int nte_types,
                                           int ni, int nj, int nx, int ny) :
    R12IntsAcc(nte_types, ni, nj, nx, ny), datafile_(MPI_FILE_NULL)
{
  filename_ = strdup(filename);

  init(false);
}

R12IntsAcc_MPIIOFile::R12IntsAcc_MPIIOFile(StateIn& si) :
  SavableState(si), R12IntsAcc(si)
{
  si.getstring(filename_);

  init(true);
}

R12IntsAcc_MPIIOFile::~R12IntsAcc_MPIIOFile() {
  for (int i = 0; i < ni(); i++)
    for (int j = 0; j < nj(); j++) {
      int ij = ij_index(i, j);
      for (int oper_type = 0; oper_type < num_te_types(); oper_type++)
        if (pairblk_[ij].ints_[oper_type] != NULL) {
          ExEnv::outn() << indent << me() << ": i = " << i << " j = "
              << j << " oper_type = " << oper_type << endl;
          throw ProgrammingError("R12IntsAcc_MPIIOFile::~R12IntsAcc_MPIIOFile -- some nonlocal blocks have not been released!",
                                 __FILE__,__LINE__);
        }
    }
  delete[] pairblk_;
  free(filename_);
}

void
R12IntsAcc_MPIIOFile::save_data_state(StateOut& so)
{
  R12IntsAcc::save_data_state(so);
  so.putstring(filename_);
}

void
R12IntsAcc_MPIIOFile::init(bool restart)
{
  int nproc = ntasks();
  int errcod;
  errcod = MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  nints_per_block_ = nxy()*num_te_types();

  pairblk_ = new struct PairBlkInfo[ni()*nj()];
  int i, j, ij;
  for(i=0,ij=0;i<ni();i++)
    for(j=0;j<nj();j++,ij++) {
      for(int type=0; type<num_te_types(); type++) {
        pairblk_[ij].ints_[type] = NULL;
        pairblk_[ij].refcount_[type] = 0;
        }
      pairblk_[ij].offset_ = (MPI_Offset)ij*blocksize();
    }

  // Try opening/creating the file
  int amode;
  if (!restart)
    amode = MPI_MODE_CREATE | MPI_MODE_WRONLY;
  else
    amode = MPI_MODE_RDONLY;

  errcod = MPI_File_open(MPI_COMM_WORLD, filename_, amode,
                         MPI_INFO_NULL, &datafile_);
  check_error_code_(errcod);
  errcod = MPI_File_close(&datafile_);
  check_error_code_(errcod);
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
R12IntsAcc_MPIIOFile::activate()
{
  R12IntsAcc::activate();
  int errcod = MPI_File_open(MPI_COMM_WORLD, filename_, MPI_MODE_RDWR, MPI_INFO_NULL, &datafile_);
  check_error_code_(errcod);
  if (classdebug() > 0)
    ExEnv::out0() << indent << "opened file=" << filename_ << " datafile=" << datafile_ << endl;
}

void
R12IntsAcc_MPIIOFile::deactivate()
{
  int errcod = MPI_File_close(&datafile_);
  check_error_code_(errcod);
  R12IntsAcc::deactivate();
  if (classdebug() > 0)
    ExEnv::out0() << indent << "closed file=" << filename_ << " datafile=" << datafile_ << endl;
}

void
R12IntsAcc_MPIIOFile::release_pair_block(int i, int j, tbint_type oper_type) const
{
  int ij = ij_index(i,j);
  struct PairBlkInfo *pb = &pairblk_[ij];
  if (pb->refcount_[oper_type] <= 0) {
    ExEnv::outn() << indent << me() << ":refcount=0: i = " << i << " j = " << j << " tbint_type = " << oper_type << endl;
    throw ProgrammingError("R12IntsAcc_MPIIOFile::release_pair_block -- refcount is already zero!",
                           __FILE__,__LINE__);
  }
  if (pb->ints_[oper_type] != NULL && pb->refcount_[oper_type] == 1) {
    delete[] pb->ints_[oper_type];
    pb->ints_[oper_type] = NULL;
  }
  pb->refcount_[oper_type] -= 1;
}

///////////////////////////////////////////////////////////////

static ClassDesc R12IntsAcc_MPIIOFile_Ind_cd(
  typeid(R12IntsAcc_MPIIOFile_Ind),"R12IntsAcc_MPIIOFile_Ind",1,"public R12IntsAcc",
  0, 0, create<R12IntsAcc_MPIIOFile_Ind>);

R12IntsAcc_MPIIOFile_Ind::R12IntsAcc_MPIIOFile_Ind(StateIn& si) :
  SavableState(si), R12IntsAcc_MPIIOFile(si)
{
}

R12IntsAcc_MPIIOFile_Ind::R12IntsAcc_MPIIOFile_Ind(const char* filename, int num_te_types,
                                                   int ni, int nj, int nx, int ny) :
  R12IntsAcc_MPIIOFile(filename,num_te_types,ni,nj,nx,ny)
{
}

R12IntsAcc_MPIIOFile_Ind::~R12IntsAcc_MPIIOFile_Ind()
{
}

void
R12IntsAcc_MPIIOFile_Ind::save_data_state(StateOut&so)
{
  R12IntsAcc_MPIIOFile::save_data_state(so);
}

#if 0
void
R12IntsAcc_MPIIOFile_Ind::store_memorygrp(Ref<MemoryGrp>& mem, int ni, const size_t blksize)
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
  else if (ni > ni()) {
    ExEnv::out0() << "R12IntsAcc_MPIIOFile_Ind::store_memorygrp(mem,ni) called with invalid argument:" << endl <<
      "ni > R12IntsAcc_MPIIOFile_Ind::ni()" << endl;
    abort();
  }
  else if (next_orbital() + ni > ni()) {
    ExEnv::out0() << "R12IntsAcc_MPIIOFile_Ind::store_memorygrp(mem,ni) called with invalid argument:" << endl <<
      "ni+next_orbital() > R12IntsAcc_MPIIOFile_Ind::ni()" << endl;
    abort();
  }
  else {
    size_t blksize_memgrp = blksize;
    if (blksize_memgrp == 0)
      blksize_memgrp = blksize();

    // Now do some extra work to figure layout of data in MemoryGrp
    // Compute global offsets to each processor's data
    int i,j,ij;
    int me = mem->me();
    int nproc = mem->n();

    // Append the data to the file
    int errcod = MPI_File_open(MPI_COMM_WORLD, filename_, MPI_MODE_CREATE | MPI_MODE_APPEND | MPI_MODE_WRONLY, MPI_INFO_NULL, &datafile_);
    check_error_code_(errcod);

    for (int i=0; i<ni; i++)
      for (int j=0; j<nj(); j++) {
	int ij = ij_index(i,j);
	int proc = ij%nproc;
        if (proc != me) continue;
	int local_ij_index = ij/nproc;
        double *data = (double *) ((size_t)mem->localdata() + blksize_memgrp*num_te_types()*local_ij_index);
        for(int te_type=0; te_type < num_te_types(); te_type++) {
          int IJ = ij_index(i+next_orbital(),j);
          int errcod = MPI_File_seek(datafile_, pairblk_[IJ].offset_+(MPI_Offset)te_type*blksize(), MPI_SEEK_SET);
          check_error_code_(errcod);
          MPI_Status status;
          errcod = MPI_File_write(datafile_, (void *)data, nxy(), MPI_DOUBLE, &status);
          check_error_code_(errcod);
          data = (double*) ((size_t) data + blksize_memgrp);
        }
      }
    // Close the file and update the i counter
    errcod = MPI_File_close(&datafile_);
    check_error_code_(errcod);
  }

  inc_next_orbital(ni);
}

void
R12IntsAcc_MPIIOFile_Ind::restore_memorygrp(Ref<MemoryGrp>& mem, int ioffset, int ni, const size_t blksize) const
{
  // mem must be the same as mem_
  if (mem_ != mem) {
    throw ProgrammingError("R12IntsAcc_MPIIOFile_Ind::restore_memorygrp() -- mem != R12IntsAcc_MPIIOFile_Ind::mem_",__FILE__,__LINE__);
  }
  if (ni != ni()) {
    throw ProgrammingError("R12IntsAcc_MPIIOFile_Ind::restore_memorygrp() -- ni != R12IntsAcc_MPIIOFile_Ind::ni()",__FILE__,__LINE__);
  }
  throw FeatureNotImplemented("R12IntsAcc_MPIIOFile_Ind::restore_memorygrp()");
}
#endif

void R12IntsAcc_MPIIOFile_Ind::store_pair_block(int i, int j,
                                                tbint_type oper_type,
                                                const double *ints)
{
  // store blocks local to this node ONLY
  assert(is_local(i,j));

  const int nproc = ntasks();
  const int ij = ij_index(i,j);
  const PairBlkInfo* pb = &pairblk_[ij];
  const int local_ij_index = ij / nproc;

  // first, seek, then write
  if (classdebug() > 0)
    ExEnv::outn() << indent << "storing block: me=" << me() << " file=" << filename_ << " i,j=" << i << "," << j << " oper_type=" << oper_type << endl;
  int errcod = MPI_File_seek(datafile_,
                             pairblk_[ij].offset_ + (MPI_Offset) oper_type * blksize(),
                             MPI_SEEK_SET);
  check_error_code_(errcod);
  MPI_Status status;
  errcod = MPI_File_write(datafile_, (void *) ints, nxy(), MPI_DOUBLE,
                          &status);
  check_error_code_(errcod);
}

const double*
R12IntsAcc_MPIIOFile_Ind::retrieve_pair_block(int i, int j, tbint_type oper_type) const
{
  int ij = ij_index(i,j);
  struct PairBlkInfo *pb = &pairblk_[ij];
  // Always first check if it's already in memory
  if (pb->ints_[oper_type] == 0) {
    MPI_Offset offset = pb->offset_ + (MPI_Offset)oper_type*blksize();
    int errcod = MPI_File_seek(datafile_, offset, MPI_SEEK_SET);
    check_error_code_(errcod);
    double *buffer = new double[nxy()];
    MPI_Status status;
    errcod = MPI_File_read(datafile_, (void *)buffer, nxy(), MPI_DOUBLE, &status);
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
