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

#include <cassert>
#include <stdexcept>
#include <stdlib.h>
#include <util/misc/consumableresources.h>
#include <util/misc/string.h>
#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <util/misc/scexception.h>
#include <math/distarray4/distarray4_mpiiofile.h>

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
  if (!this->active()) {
    DistArray4::activate();
    int errcod = MPI_File_open(MPI_COMM_WORLD, filename_, MPI_MODE_RDWR, MPI_INFO_NULL, &datafile_);
    check_error_code_(errcod);
    if (classdebug() > 0)
      ExEnv::out0() << indent << "opened file=" << filename_ << " datafile=" << datafile_ << endl;
  }
}

void
DistArray4_MPIIOFile::deactivate()
{
  if (this->active()) {
    int errcod = MPI_File_close(&datafile_);
    check_error_code_(errcod);
    DistArray4::deactivate();
    if (classdebug() > 0)
      ExEnv::out0() << indent << "closed file=" << filename_ << " datafile=" << datafile_ << endl;
  }
}

void
DistArray4_MPIIOFile::release_pair_block(int i, int j, tbint_type oper_type) const
{
  MPQC_ASSERT(this->active());  //make sure we are active
  int ij = ij_index(i,j);
  struct PairBlkInfo *pb = &pairblk_[ij];
  if (pb->refcount_[oper_type] <= 0) {
    ExEnv::outn() << indent << me() << ":refcount=0: i = " << i << " j = " << j << " tbint_type = " << oper_type << endl;
    throw ProgrammingError("DistArray4_MPIIOFile::release_pair_block -- refcount is already zero!",
                           __FILE__,__LINE__);
  }
  if (pb->ints_[oper_type] != NULL && pb->refcount_[oper_type] == 1) {
    if (pb->manage_[oper_type]) { // deallocate only if I am in charge
      deallocate(pb->ints_[oper_type]);
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

#ifndef AVOID_XLC_BUG
Ref<DistArray4>
DistArray4_MPIIOFile_Ind::clone(const DistArray4Dimensions& dim) {
  return DistArray4_MPIIOFile::clone<DistArray4_MPIIOFile_Ind>(dim);
}
#else
Ref<DistArray4>
DistArray4_MPIIOFile_Ind::clone(const DistArray4Dimensions& dim) {

    int id = 0;
    std::string clonename;
    using detail::clone_filename;
    clone_filename(clonename, this->filename_, id);
    if (clonelist_.nonnull()) {
      while (clonelist_->key_exists(clonename)) {
        ++id;
        clone_filename(clonename, this->filename_, id);
      }
    } else {
      clonelist_ = ListOfClones::instance();
    }
    clonelist_->add(clonename, id);

    Ref<DistArray4_MPIIOFile_Ind> result;
    if (dim == DistArray4Dimensions::default_dim())
      result = new DistArray4_MPIIOFile_Ind(clonename.c_str(), num_te_types(),
                           ni(), nj(),
                           nx(), ny(),
                           storage());
    else
      result = new DistArray4_MPIIOFile_Ind(clonename.c_str(), dim.num_te_types(),
                           dim.n1(), dim.n2(),
                           dim.n3(), dim.n4(),
                           dim.storage());

    result->set_clonelist(clonelist_);
    return result;
  }
#endif

void DistArray4_MPIIOFile_Ind::store_pair_block(int i, int j,
                                                tbint_type oper_type,
                                                const double *ints)
{
  MPQC_ASSERT(this->active());  //make sure we are active

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

void DistArray4_MPIIOFile_Ind::store_pair_subblock(int i, int j, tbint_type oper_type,
                                                   int xstart, int xfence, int ystart, int yfence,
                                                   const double *buf)
{
  MPQC_ASSERT(this->active());  //make sure we are active

  const bool contiguous = (ystart == 0) && (yfence == ny());
  const int xsize = xfence - xstart;
  const int ysize = yfence - ystart;
  const int xysize = xsize * ysize;
  const int bufsize = xysize * sizeof(double);

  const int nproc = ntasks();
  const int ij = ij_index(i,j);
  const PairBlkInfo* pb = &pairblk_[ij];
  const int local_ij_index = ij / nproc;

  const size_t batchsize = contiguous ? xysize : ysize;
  const MPI_Offset stridesize = (MPI_Offset)ny() * sizeof(double);
  MPI_Offset offset = pairblk_[ij].offset_ + (MPI_Offset) oper_type * blksize() +
                                             (MPI_Offset) (xstart*ny() + ystart)*sizeof(double);
  size_t wrote_this_many = 0;

  while (wrote_this_many < xysize) {
    // first, seek, then write
    if (classdebug() > 0)
      ExEnv::outn() << indent << "storing block: me=" << me() << " file=" << filename_ << " i,j=" << i << "," << j
                    << " oper_type=" << oper_type << " batchsize=" << batchsize << endl;
    int errcod = MPI_File_seek(datafile_,
                               offset,
                               MPI_SEEK_SET);
    check_error_code_(errcod);
    MPI_Status status;
    errcod = MPI_File_write(datafile_, (void *) buf, batchsize, MPI_DOUBLE,
                            &status);
    check_error_code_(errcod);
    wrote_this_many += batchsize;
    offset += stridesize;
    buf += batchsize;
  }
}

const double*
DistArray4_MPIIOFile_Ind::retrieve_pair_block(int i, int j, tbint_type oper_type,
                                              double* buf) const
{
  MPQC_ASSERT(this->active());  //make sure we are active
  int ij = ij_index(i,j);
  struct PairBlkInfo *pb = &pairblk_[ij];
  // Always first check if it's already in memory
  // if this block is not managed, assume that user takes care of the details
  // hence read it in again
  if (pb->ints_[oper_type] == 0 || pb->manage_[oper_type] == false) {
    MPI_Offset offset = pb->offset_ + (MPI_Offset)oper_type*blksize();
    int errcod = MPI_File_seek(datafile_, offset, MPI_SEEK_SET);
    check_error_code_(errcod);
    double *buffer;
    if (buf == 0) {
      buffer = allocate<double>(nxy());
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

namespace {
  template <typename T>
  class ScratchBuffer {
    public:
    ScratchBuffer() : size_(0), buf_(0) {}
    ~ScratchBuffer() { if (buf_) free(buf_); buf_ = 0; size_ = 0;}
    T* buffer(size_t size) {
      if (size <= size_) return buf_;
      buf_ = (T*)realloc((void*)buf_,size);
      size_ = size;
      return buf_;
    }

    private:
      size_t size_;
      T* buf_;
  };
}

void
DistArray4_MPIIOFile_Ind::retrieve_pair_subblock(int i, int j, tbint_type oper_type,
                                                 int xstart, int xfence, int ystart, int yfence,
                                                 double* buf) const
{
  MPQC_ASSERT(this->active());  //make sure we are active
  static ScratchBuffer<char> scratch;
  Ref<ThreadLock> scratch_lock = ThreadGrp::get_default_threadgrp()->new_lock();

  const bool contiguous = (ystart == 0) && (yfence == ny());
  const int xsize = xfence - xstart;
  const int ysize = yfence - ystart;
  const int xysize = xsize * ysize;
  const int bufsize = xysize * sizeof(double);

  int ij = ij_index(i,j);
  struct PairBlkInfo *pb = &pairblk_[ij];
  // Always first check if it's already in memory
  // if this block is not managed, assume that user takes care of the details
  // hence read it in again
  if (pb->ints_[oper_type] == 0 || pb->manage_[oper_type] == false) {
    MPI_Offset offset = pb->offset_ + (MPI_Offset)oper_type*blksize() +
                        (MPI_Offset)(xstart*ny() + ystart)*sizeof(double);
    int errcod = MPI_File_seek(datafile_, offset, MPI_SEEK_SET);
    check_error_code_(errcod);

    // do not assume that there is enough memory -- use static scratch and lock, if necessary
    size_t readbuf_size = contiguous ? bufsize : ((xsize-1) * ny() + ysize) * sizeof(double);
    void* readbuf = buf;
    if (!contiguous) {
      scratch_lock->lock();
      readbuf = scratch.buffer(readbuf_size);
    }

    // read
    MPI_Status status;
    errcod = MPI_File_read(datafile_, (void *)readbuf, readbuf_size, MPI_DOUBLE, &status);
    check_error_code_(errcod);

    // hence may need to copy the required data to the buffer
    if (!contiguous) {
      double* outbuf = buf;
      const double* srcbuf = (const double*) readbuf;
      for(int x=0; x<xsize; ++x, srcbuf+=ny(), outbuf+=ysize) {
        std::copy(srcbuf, srcbuf + ysize, outbuf);
      }
      // release scratch if necessary
      if (!contiguous) scratch_lock->unlock();
    }
  }
  else { // if data is already available need to copy to the buffer
    double* outbuf = buf;
    const double* srcbuf = pb->ints_[oper_type] + (xstart * ny() + ystart);
    if (contiguous) { // read all at once
      std::copy(srcbuf, srcbuf + xysize, outbuf);
    }
    else { // read row by row
      for(int x=0; x<xsize; ++x, srcbuf+=ny(), outbuf+=ysize) {
        std::copy(srcbuf, srcbuf + ysize, outbuf);
      }
    }
  }
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
