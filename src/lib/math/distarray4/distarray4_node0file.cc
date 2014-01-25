//
// distarray4_node0file.cc
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
#include <util/misc/consumableresources.h>
#include <util/misc/scexception.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <stdlib.h>
#include <util/misc/string.h>
#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <math/distarray4/distarray4_node0file.h>

#define CREATE_FILE_ON_NODE0_ONLY 1

using namespace std;
using namespace sc;

///////////////////////////////////////////////////////////////

static ClassDesc DistArray4_Node0File_cd(
  typeid(DistArray4_Node0File),"DistArray4_Node0File",1,"public DistArray4",
  0, 0, create<DistArray4_Node0File>);

DistArray4_Node0File::DistArray4_Node0File(const char* filename, int num_te_types,
                                           int ni, int nj, int nx, int ny,
                                           DistArray4Storage storage) :
  DistArray4(num_te_types, ni, nj, nx, ny, storage)
{
  filename_ = strdup(filename);
  init(false);
}

DistArray4_Node0File::DistArray4_Node0File(StateIn& si) :
  DistArray4(si)
{
  si.getstring(filename_);
  clonelist_ = ListOfClones::restore_instance(si);
  init(true);
}

DistArray4_Node0File::~DistArray4_Node0File() {
  for (int i = 0; i < ni(); i++)
    for (int j = 0; j < nj(); j++) {
      if (!is_avail(i, j)) {
        int ij = ij_index(i, j);
        for (int oper_type = 0; oper_type < num_te_types(); oper_type++)
          if (pairblk_[ij].ints_[oper_type] != NULL) {
            ExEnv::outn() << indent << me() << ": i = " << i << " j = " << j
                << " oper_type = " << oper_type << endl;
            throw ProgrammingError("DistArray4_Node0File::~DistArray4_Node0File -- some nonlocal blocks have not been released!",
                                   __FILE__,__LINE__);
          }
      }
    }
  delete[] pairblk_;

  // Destroy the file
#if  CREATE_FILE_ON_NODE0_ONLY
  if (me() == 0)
#endif
    unlink(filename_);
  free(filename_);
}

void
DistArray4_Node0File::save_data_state(StateOut& so)
{
  DistArray4::save_data_state(so);
  so.putstring(filename_);
  ListOfClones::save_instance(clonelist_, so);
}

namespace {
  void clone_filename(std::string& result, const char* original, int id) {
    std::ostringstream oss;
    oss << original << ".clone" << id;
    result = oss.str();
  }
}

Ref<DistArray4>
DistArray4_Node0File::clone(const DistArray4Dimensions& dim) {
  int id = 0;
  std::string clonename;
  clone_filename(clonename, this->filename_, id);
  if (clonelist_.nonnull()) {
    while (clonelist_->key_exists(clonename)) {
      ++id;
      clone_filename(clonename, this->filename_, id);
    }
  }
  else {
    clonelist_ = ListOfClones::instance();
  }
  clonelist_->add(clonename, id);

  Ref<DistArray4_Node0File> result;
  if (dim == DistArray4Dimensions::default_dim())
    result = new DistArray4_Node0File(clonename.c_str(), num_te_types(),
                                      ni(), nj(), nx(), ny(), storage());
  else
    result = new DistArray4_Node0File(clonename.c_str(), dim.num_te_types(),
                                      dim.n1(), dim.n2(), dim.n3(), dim.n4(),
                                      dim.storage());

  result->set_clonelist(clonelist_);
  return result;
}

void
DistArray4_Node0File::init(bool restart)
{
  pairblk_ = new PairBlkInfo[ni()*nj()];
  int i, j, ij;
  for(i=0,ij=0;i<ni();i++)
    for(j=0;j<nj();j++,ij++) {
      for(int type=0; type<num_te_types(); type++) {
        pairblk_[ij].ints_[type] = NULL;
        pairblk_[ij].manage_[type] = false;
        pairblk_[ij].refcount_[type] = 0;
        if (classdebug() > 0)
          ExEnv::outn() << indent << me() << ":refcount=" << pairblk_[ij].refcount_[type]
                        << ": i = " << i << " j = " << j << " tbint_type = " << type << endl;
        }
      pairblk_[ij].offset_ = (off_t)ij*blocksize();
    }

  // node 0 will have the file
#if CREATE_FILE_ON_NODE0_ONLY
  if (me() != 0)
    return;
#endif

  // See if can open/create the file
  if (restart) {
    datafile_ = open(filename_,O_RDWR,0644);
  }
  else {
    datafile_ = creat(filename_,0644);
  }
  // Check if the file was opened correctly
  check_filedescr_();
  // If everything is fine close it and proceed
  close(datafile_);
}

void
DistArray4_Node0File::check_filedescr_()
{
  // If the file was not opened correctly - throw an exception
  if (datafile_ == -1) {
    switch (errno) {
      case EACCES:
      throw std::runtime_error("DistArray4_Node0File::DistArray4_Node0File -- access to the requested file is not allowed");
      case ENOSPC:
      throw std::runtime_error("DistArray4_Node0File::DistArray4_Node0File -- no space left in the filesystem");
      case ENOENT:
      throw std::runtime_error("DistArray4_Node0File::DistArray4_Node0File -- file does not exist");
      case ENOTDIR:
      throw std::runtime_error("DistArray4_Node0File::DistArray4_Node0File -- not a directory");
      default:
      const char* errormsg = strerror(errno);
      ExEnv::out0() << "DistArray4_Node0File::DistArray4_Node0File " << errormsg << std::endl;
      throw std::runtime_error("DistArray4_Node0File::DistArray4_Node0File -- failed to open POSIX file on node 0");
    }
  }
}

void
DistArray4_Node0File::set_clonelist(const Ref<ListOfClones>& cl) {
  clonelist_ = cl;
}

void
DistArray4_Node0File::activate()
{
  if (active()) return;

#if CREATE_FILE_ON_NODE0_ONLY
  if (me() == 0)
#endif
    datafile_ = open(filename_, O_RDWR);
  DistArray4::activate();
  if (classdebug() > 0)
    ExEnv::out0() << indent << "opened file=" << filename_ << " datafile=" << datafile_ << endl;
}

void
DistArray4_Node0File::deactivate()
{
  if (!active()) return;

#if CREATE_FILE_ON_NODE0_ONLY
  if (me() == 0)
#endif
    close(datafile_);
  DistArray4::deactivate();
  if (classdebug() > 0)
    ExEnv::out0() << indent << "closed file=" << filename_ << " datafile=" << datafile_ << endl;
}

void
DistArray4_Node0File::store_pair_block(int i, int j, tbint_type oper_type, const double *data)
{
  MPQC_ASSERT(this->active());  //make sure we are active
  // Can write blocks?
  if (!is_avail(i,j))
    throw ProgrammingError("DistArray4_Node0File::store_pair_block -- can only be called on node 0",
                           __FILE__,__LINE__);

  const int ij = ij_index(i,j);
  const PairBlkInfo* pb = &pairblk_[ij];

  // first, seek
  const off_t offset = pb->offset_ + (off_t)oper_type*blksize();
  if (classdebug() > 0)
    ExEnv::out0() << indent << "storing block: file=" << filename_ << " i,j=" << i << "," << j << " oper_type=" << oper_type << " offset=" << offset << endl;
  const off_t result_offset = lseek(datafile_,offset,SEEK_SET);
  if (offset == static_cast<off_t>(-1) || result_offset != offset) {
    std::ostringstream oss;
    oss << "DistArray4_Node0File::store_pair_block() -- lseek failed: " << strerror(errno);
    throw FileOperationFailed(oss.str().c_str(),
                                __FILE__, __LINE__,
                                filename_, FileOperationFailed::Other);
  }

  // then, write
  ssize_t wrote_this_much = write(datafile_, data, blksize());
  if (wrote_this_much != blksize()) {
    const char* errormsg = strerror(errno);
    ExEnv::out0() << "DistArray4_Node0File::store_pair_block(): " << errormsg << std::endl;
    throw FileOperationFailed("DistArray4_Node0File::store_pair_block() -- write failed",
                              __FILE__, __LINE__,
                              filename_, FileOperationFailed::Write);
  }
}

void
DistArray4_Node0File::store_pair_subblock(int i, int j, tbint_type oper_type,
                                          int xstart, int xfence, int ystart, int yfence,
                                          const double *buf)
{
  MPQC_ASSERT(this->active());  //make sure we are active
  // Can write blocks?
  if (!is_avail(i,j))
    throw ProgrammingError("DistArray4_Node0File::store_pair_block -- can only be called on node 0",
                           __FILE__,__LINE__);

  const bool contiguous = (ystart == 0) && (yfence == ny());
  const int xsize = xfence - xstart;
  const int ysize = yfence - ystart;
  const int xysize = xsize * ysize;
  const int bufsize = xysize * sizeof(double);

  const int ij = ij_index(i,j);
  const PairBlkInfo* pb = &pairblk_[ij];

  const size_t batchsize = (contiguous ? xysize : ysize) * sizeof(double);
  const off_t stridesize = (off_t)ny() * sizeof(double);
  off_t offset = pb->offset_ +
          (off_t)oper_type*blksize() +
          (off_t)(xstart*ny() + ystart)*sizeof(double);
  ssize_t wrote_this_much = 0;
  while (wrote_this_much < bufsize) {
    // first, seek
    if (classdebug() > 0)
      ExEnv::out0() << indent << "storing block: file=" << filename_ << " i,j=" << i << "," << j << " oper_type=" << oper_type
                    << " offset=" << offset << " batchsize=" << batchsize << endl;
    const off_t result_offset = lseek(datafile_,offset,SEEK_SET);
    if (offset == static_cast<off_t>(-1) || result_offset != offset) {
      std::ostringstream oss;
      oss << "DistArray4_Node0File::store_pair_block() -- lseek failed: " << strerror(errno);
      throw FileOperationFailed(oss.str().c_str(),
                                __FILE__, __LINE__,
                                filename_, FileOperationFailed::Other);
    }

    // then, write
    const ssize_t amount = write(datafile_, buf, batchsize);
    wrote_this_much += amount;
    if (amount != batchsize) {
      const char* errormsg = strerror(errno);
      ExEnv::out0() << "DistArray4_Node0File::store_pair_block(): " << errormsg << std::endl;
      throw FileOperationFailed("DistArray4_Node0File::store_pair_block() -- write failed",
                                __FILE__, __LINE__,
                                filename_, FileOperationFailed::Write);
    }
    offset += stridesize;
    buf += batchsize/sizeof(double);
  }
}

const double * DistArray4_Node0File::retrieve_pair_block(int i, int j,
                                                   tbint_type oper_type,
                                                   double* buf) const {
  if (not this->active()) { //make sure we are active
    std::ostringstream oss;
    oss << "DistArray4_Node0File::retrieve_pair_block -- file " << this->filename_ << " is not open" << std::endl;
    ExEnv::outn() << oss.str();
    throw ProgrammingError(oss.str().c_str(), __FILE__, __LINE__);
  }
  // Can read blocks?
  if (!is_avail(i, j))
    throw ProgrammingError("DistArray4_Node0File::retrieve_pair_block -- can only be called on node 0",
        __FILE__,__LINE__);

  const int ij = ij_index(i, j);
  const PairBlkInfo* pb = &pairblk_[ij];
  // Always first check if it's already in memory
  // if I don't manage the memory for this block, assume that the user is in charge of memory management
  // therefore read in again
  if (pb->ints_[oper_type] == 0 || pb->manage_[oper_type] == false) {

    if (classdebug() > 0)
      ExEnv::out0() << indent << "retrieving block: file=" << filename_
          << " i,j=" << i << "," << j << " oper_type=" << oper_type << endl;

    off_t offset = pb->offset_ + (off_t)oper_type*blksize();
    off_t result_offset = lseek(datafile_, offset, SEEK_SET);
    if (offset == (off_t)-1 || result_offset != offset) {
      std::ostringstream oss;
      oss << "DistArray4_Node0File::store_pair_block() -- lseek failed: " << strerror(errno);
      throw FileOperationFailed(oss.str().c_str(),
          __FILE__,
          __LINE__,
          filename_,
          FileOperationFailed::Other);
    }
    if (buf != 0) {
      pb->ints_[oper_type] = buf;
      pb->manage_[oper_type] = false;
    }
    else {
      pb->ints_[oper_type] = allocate<double>(nxy());
      pb->manage_[oper_type] = true;
    }
    ssize_t read_this_much = read(datafile_, pb->ints_[oper_type], blksize());
    if (read_this_much != blksize()) {
      std::ostringstream oss;
      oss << "DistArray4_Node0File::store_pair_block() -- read failed: " << strerror(errno);
      throw FileOperationFailed(oss.str().c_str(),
          __FILE__,
          __LINE__,
          filename_,
          FileOperationFailed::Read);
    }
  }
  else { // data is already available
    if (buf != 0 && buf != pb->ints_[oper_type]) // may need to copy
      std::copy(pb->ints_[oper_type], pb->ints_[oper_type] + this->nxy(), buf);
  }
  pb->refcount_[oper_type] += 1;
  if (classdebug() > 0)
    ExEnv::outn() << indent << me() << ":refcount="
        << pb->refcount_[oper_type] << ": i = " << i << " j = " << j
        << " tbint_type = " << oper_type << endl;
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
DistArray4_Node0File::retrieve_pair_subblock(int i, int j, tbint_type oper_type,
                                             int xstart, int xfence, int ystart, int yfence,
                                             double* buf) const
{
  MPQC_ASSERT(this->active());  //make sure we are active
  static ScratchBuffer<char> scratch;
  static Ref<ThreadLock> read_lock = ThreadGrp::get_default_threadgrp()->new_lock();

  const int xsize = xfence - xstart;
  const int ysize = yfence - ystart;
  const int xysize = xsize * ysize;
  const int bufsize = xysize * sizeof(double);
  const bool contiguous = (ysize == ny()) || (xsize == 1);

  // Can read blocks?
  if (!is_avail(i, j))
    throw ProgrammingError("DistArray4_Node0File::retrieve_pair_block() -- can only be called on node 0",
        __FILE__,__LINE__);

  const int ij = ij_index(i, j);
  const PairBlkInfo* pb = &pairblk_[ij];
  // Always first check if it's already in memory
  // if I don't manage the memory for this block, assume that the user is in charge of memory management
  // therefore need to read in again
  if (pb->ints_[oper_type] == 0 || pb->manage_[oper_type] == false) {

    off_t offset = pb->offset_ + (off_t)oper_type*blksize() +
                   (off_t)(xstart*ny() + ystart)*sizeof(double);
    read_lock->lock();
    off_t result_offset = lseek(datafile_, offset, SEEK_SET);
    if (offset == (off_t)-1 || result_offset != offset) {
      std::ostringstream oss;
      oss << "DistArray4_Node0File::retrieve_pair_block() -- lseek failed: " << strerror(errno);
      read_lock->unlock();
      throw FileOperationFailed(oss.str().c_str(),
          __FILE__,
          __LINE__,
          filename_,
          FileOperationFailed::Other);
    }

    // do not assume that there is enough memory -- use static scratch and lock, if necessary
    size_t readbuf_size = contiguous ? bufsize : ((xsize-1) * ny() + ysize) * sizeof(double);
    void* readbuf = buf;
    if (!contiguous) {
      readbuf = scratch.buffer(readbuf_size);
    }

    const ssize_t read_this_much = read(datafile_, readbuf, readbuf_size);
    if (read_this_much != readbuf_size) {
      std::ostringstream oss;
      oss << "DistArray4_Node0File::retrieve_pair_subblock() -- read failed: " << strerror(errno);
      read_lock->unlock();
      throw FileOperationFailed(oss.str().c_str(),
                                __FILE__,
                                __LINE__,
                                filename_,
                                FileOperationFailed::Read);
    }

    // hence may need to copy the required data to the buffer
    if (!contiguous) {
      double* outbuf = buf;
      const double* srcbuf = (const double*) readbuf;
      for(int x=0; x<xsize; ++x, srcbuf+=ny(), outbuf+=ysize) {
        std::copy(srcbuf, srcbuf + ysize, outbuf);
      }
    }
    read_lock->unlock();
  }
  else { // data is already available
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

void
DistArray4_Node0File::release_pair_block(int i, int j, tbint_type oper_type) const
{
  MPQC_ASSERT(this->active());  //make sure we are active
  if (is_avail(i,j)) {
    const int ij = ij_index(i,j);
    const PairBlkInfo *pb = &pairblk_[ij];
    if (pb->refcount_[oper_type] <= 0) {
      ExEnv::outn() << indent << me() << ":refcount=0: i = " << i << " j = " << j << " tbint_type = " << oper_type << endl;
      throw std::runtime_error("Logic error: DistArray4_Node0File::release_pair_block: refcount is already zero!");
    }
    if (pb->ints_[oper_type] != NULL && pb->refcount_[oper_type] == 1) {
      if (classdebug() > 0)
        ExEnv::out0() << indent << "releasing block: file=" << filename_ << " i,j=" << i << "," << j << " oper_type=" << oper_type << endl;
      if (pb->manage_[oper_type]) // deallocate if managed by me
        deallocate(pb->ints_[oper_type]);
      pb->ints_[oper_type] = NULL;
      pb->manage_[oper_type] = false;
    }
    pb->refcount_[oper_type] -= 1;
    if (classdebug() > 0)
      ExEnv::outn() << indent << me() << ":refcount=" << pb->refcount_[oper_type]
                    << ": i = " << i << " j = " << j << " tbint_type = " << oper_type << endl;
  }
}


// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
