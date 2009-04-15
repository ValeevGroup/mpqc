//
// distarray4_mpiiofile.cc
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
#include <chemistry/qc/mbptr12/distarray4_mpiiofile.h>

using namespace std;
using namespace sc;

///////////////////////////////////////////////////////////////

static ClassDesc DistArray4_MPIIOFile_cd(
  typeid(DistArray4_MPIIOFile),"DistArray4_MPIIOFile",1,"public DistArray4",
  0, 0, 0);

DistArray4_MPIIOFile::DistArray4_MPIIOFile(const char* filename, int nte_types,
                                           int ni, int nj, int nx, int ny,
                                           DistArray4Storage storage) :
    DistArray4(nte_types, ni, nj, nx, ny, storage), datafile_(MPI_FILE_NULL)
{
  filename_ = strdup(filename);

  init(false);
}

DistArray4_MPIIOFile::DistArray4_MPIIOFile(StateIn& si) :
  SavableState(si), DistArray4(si)
{
  si.getstring(filename_);
  clonelist_ = ListOfClones::restore_instance(si);

  init(true);
}

DistArray4_MPIIOFile::~DistArray4_MPIIOFile() {
  for (int i = 0; i < ni(); i++)
    for (int j = 0; j < nj(); j++) {
      int ij = ij_index(i, j);
      for (int oper_type = 0; oper_type < num_te_types(); oper_type++)
        if (pairblk_[ij].ints_[oper_type] != NULL) {
          ExEnv::outn() << indent << me() << ": i = " << i << " j = "
              << j << " oper_type = " << oper_type << endl;
          throw ProgrammingError("DistArray4_MPIIOFile::~DistArray4_MPIIOFile -- some nonlocal blocks have not been released!",
                                 __FILE__,__LINE__);
        }
    }
  delete[] pairblk_;
  free(filename_);
}

void
DistArray4_MPIIOFile::save_data_state(StateOut& so)
{
  DistArray4::save_data_state(so);
  so.putstring(filename_);
  ListOfClones::save_instance(clonelist_,so);
}

void
DistArray4_MPIIOFile::init(bool restart)
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
        pairblk_[ij].manage_[type] = false;
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
DistArray4_MPIIOFile::check_error_code_(int errcod) const
{

  if (errcod != MPI_SUCCESS) {

    int errclass;
    MPI_Error_class(errcod, &errclass);

    switch (errclass) {
    default:
      char* errstr = new char[MPI_MAX_ERROR_STRING];
      int errstrlen;
      MPI_Error_string(errcod, errstr, &errstrlen);
      ExEnv::out0() << "DistArray4_MPIIOFile::DistArray4_MPIIOFile -- MPI-I/O error: " << errstr
                    << " on file " << filename_ << endl;
      delete[] errstr;
      throw std::runtime_error("DistArray4_MPIIOFile::DistArray4_MPIIOFile -- MPI-I/O error");
    }
  }

}

void
DistArray4_MPIIOFile::set_clonelist(const Ref<ListOfClones>& cl) {
  clonelist_ = cl;
}

void
DistArray4_MPIIOFile::activate()
{
  DistArray4::activate();
  int errcod = MPI_File_open(MPI_COMM_WORLD, filename_, MPI_MODE_RDWR, MPI_INFO_NULL, &datafile_);
  check_error_code_(errcod);
  if (classdebug() > 0)
    ExEnv::out0() << indent << "opened file=" << filename_ << " datafile=" << datafile_ << endl;
}

void
DistArray4_MPIIOFile::deactivate()
{
  int errcod = MPI_File_close(&datafile_);
  check_error_code_(errcod);
  DistArray4::deactivate();
  if (classdebug() > 0)
    ExEnv::out0() << indent << "closed file=" << filename_ << " datafile=" << datafile_ << endl;
}

void
DistArray4_MPIIOFile::release_pair_block(int i, int j, tbint_type oper_type) const
{
  int ij = ij_index(i,j);
  struct PairBlkInfo *pb = &pairblk_[ij];
  if (pb->refcount_[oper_type] <= 0) {
    ExEnv::outn() << indent << me() << ":refcount=0: i = " << i << " j = " << j << " tbint_type = " << oper_type << endl;
    throw ProgrammingError("DistArray4_MPIIOFile::release_pair_block -- refcount is already zero!",
                           __FILE__,__LINE__);
  }
  if (pb->ints_[oper_type] != NULL && pb->refcount_[oper_type] == 1) {
    if (pb->manage_[oper_type]) { // deallocate only if I am in charge
      delete[] pb->ints_[oper_type];
      pb->manage_[oper_type] = false;
    }
    pb->ints_[oper_type] = NULL;
  }
  pb->refcount_[oper_type] -= 1;
}

///////////////////////////////////////////////////////////////

static ClassDesc DistArray4_MPIIOFile_Ind_cd(
  typeid(DistArray4_MPIIOFile_Ind),"DistArray4_MPIIOFile_Ind",1,"public DistArray4",
  0, 0, create<DistArray4_MPIIOFile_Ind>);

DistArray4_MPIIOFile_Ind::DistArray4_MPIIOFile_Ind(StateIn& si) :
  SavableState(si), DistArray4_MPIIOFile(si)
{
}

DistArray4_MPIIOFile_Ind::DistArray4_MPIIOFile_Ind(const char* filename, int num_te_types,
                                                   int ni, int nj, int nx, int ny,
                                                   DistArray4Storage storage) :
  DistArray4_MPIIOFile(filename,num_te_types,ni,nj,nx,ny, storage)
{
}

DistArray4_MPIIOFile_Ind::~DistArray4_MPIIOFile_Ind()
{
}

void
DistArray4_MPIIOFile_Ind::save_data_state(StateOut&so)
{
  DistArray4_MPIIOFile::save_data_state(so);
}

Ref<DistArray4>
DistArray4_MPIIOFile_Ind::clone(const DistArray4Dimensions& dim) {
  return DistArray4_MPIIOFile::clone<DistArray4_MPIIOFile_Ind>(dim);
}

void DistArray4_MPIIOFile_Ind::store_pair_block(int i, int j,
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
DistArray4_MPIIOFile_Ind::retrieve_pair_block(int i, int j, tbint_type oper_type,
                                              double* buf) const
{
  int ij = ij_index(i,j);
  struct PairBlkInfo *pb = &pairblk_[ij];
  // Always first check if it's already in memory
  if (pb->ints_[oper_type] == 0) {
    MPI_Offset offset = pb->offset_ + (MPI_Offset)oper_type*blksize();
    int errcod = MPI_File_seek(datafile_, offset, MPI_SEEK_SET);
    check_error_code_(errcod);
    double *buffer;
    if (buf == 0) {
      buffer = new double[nxy()];
      pb->manage_[oper_type] = true;
    }
    else {
      buffer = buf;
      pb->manage_[oper_type] = false;
    }
    MPI_Status status;
    errcod = MPI_File_read(datafile_, (void *)buffer, nxy(), MPI_DOUBLE, &status);
    check_error_code_(errcod);
    pb->ints_[oper_type] = buffer;
  }
  else { // if data is already available may need to copy the buffer
    if (buf != 0 && buf != pb->ints_[oper_type])
      std::copy(pb->ints_[oper_type], pb->ints_[oper_type]+this->nxy(), buf);
  }
  pb->refcount_[oper_type] += 1;
  if (buf)
    return buf;
  else
    return pb->ints_[oper_type];
}

void detail::clone_filename(std::string& result, const char* original, int id) {
  std::ostringstream oss;
  oss << original << ".clone" << id;
  result = oss.str();
}

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
