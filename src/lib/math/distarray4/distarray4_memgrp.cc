//
// distarray4_memgrp.cc
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

#include <stdexcept>
#include <cassert>
#include <cstdlib>
#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <util/misc/scexception.h>
#include <math/distarray4/distarray4_memgrp.h>

using namespace std;
using namespace sc;

///////////////////////////////////////////////////////////////

static ClassDesc DistArray4_MemoryGrp_cd(
  typeid(DistArray4_MemoryGrp),"DistArray4_MemoryGrp",1,"public DistArray4",
  0, 0, create<DistArray4_MemoryGrp>);

DistArray4_MemoryGrp::DistArray4_MemoryGrp(const Ref<MemoryGrp>& mem, int num_te_types,
                                           int ni, int nj, int nx, int ny,
                                           size_t blksize_memgrp,
                                           DistArray4Storage storage) :
  DistArray4(num_te_types, ni, nj, nx, ny, storage), blksize_memgrp_(blksize_memgrp), mem_(mem)
{
  init();
}

DistArray4_MemoryGrp::DistArray4_MemoryGrp(StateIn& si) : DistArray4(si), mem_(MemoryGrp::get_default_memorygrp())
{
  si.get(blksize_memgrp_);

  init();
}

DistArray4_MemoryGrp::~DistArray4_MemoryGrp() {
  for (int i=0; i<ni(); i++)
    for (int j=0; j<nj(); j++)
      if (!is_local(i, j)) {
        int ij = ij_index(i, j);
        for (int oper_type=0; oper_type<num_te_types(); oper_type++)
          if (pairblk_[ij].ints_[oper_type] != NULL) {
            ExEnv::outn() << indent << me() << ": i = " << i << " j = "
                << j << " oper_type = " << oper_type << endl;
            throw ProgrammingError("Logic error: DistArray4_MemoryGrp::~ : some nonlocal blocks have not been released!",
                                   __FILE__, __LINE__);
          }
      }
  delete[] pairblk_;
}

void
DistArray4_MemoryGrp::save_data_state(StateOut& so)
{
  DistArray4::save_data_state(so);
  so.put(blksize_memgrp_);
}

Ref<DistArray4>
DistArray4_MemoryGrp::clone(const DistArray4Dimensions& dim)
{
  Ref<MemoryGrp> newmem = mem_->clone();
  newmem->set_localsize(mem_->localsize());
  Ref<DistArray4> result;
  if (dim == DistArray4Dimensions::default_dim())
    result = new DistArray4_MemoryGrp(newmem, num_te_types(),
                                      ni(), nj(), nx(), ny(),
                                      blksize_memgrp_,
                                      storage());
  else
    result = new DistArray4_MemoryGrp(newmem, dim.num_te_types(),
                                      dim.n1(), dim.n2(), dim.n3(), dim.n4(),
                                      dim.n3() * dim.n4() * sizeof(double),
                                      dim.storage());

  return result;
}

void DistArray4_MemoryGrp::init() {
  const int n = ntasks();

  // Now do some extra work to figure layout of data in MemoryGrp
  // Compute global offsets to each processor's data
  int i, j, ij;
  pairblk_ = new struct PairBlkInfo[ni()*nj()];
  for (i=0, ij=0; i<ni(); i++)
    for (j=0; j<nj(); j++, ij++) {
      if (not is_local(i,j)) { // remote blocks are not yet cached here
        for (int type=0; type<num_te_types(); type++) {
          pairblk_[ij].ints_[type] = NULL;
          pairblk_[ij].refcount_[type] = 0;
        }
      }
      else { // location of local blocks is immediately known
        const int ij = ij_index(i,j);
        const int local_ij_index = ij / ntasks();
        for (int type=0; type<num_te_types(); type++) {
          pairblk_[ij].ints_[type] =
              const_cast<const double *>(reinterpret_cast<double*>((size_t)mem_->localdata() +
                  blksize_memgrp_*(num_te_types()*local_ij_index + type)));
          pairblk_[ij].refcount_[type] = 0;
        }
      }
    }
}

void
DistArray4_MemoryGrp::store_pair_block(int i, int j, tbint_type oper_type, const double *ints)
{
  MPQC_ASSERT(this->active());  //make sure we are active
  MPQC_ASSERT(i>=0 && i<ni());
  MPQC_ASSERT(j>=0 && j<nj());

  const int ij = ij_index(i,j);
  struct PairBlkInfo *pb = &pairblk_[ij];
  // whether it's local or not writing will happen via the MemoryGrp
  const int local_ij_index = ij / ntasks();
  pb->offset_ = mem_->offset(ij_proc(i,j)) + (distsize_t)local_ij_index*blksize_memgrp_*num_te_types();
  const distsize_t offset = pb->offset_ + (distsize_t)oper_type*blksize_memgrp_;
  mem_->write(static_cast<const void*>(ints), offset, this->blksize());
  // cannot copy directly since may need to lock ... MemoryGrp may need a write() method.
  //void* dest = mem_->obtain_writeonly(offset, blksize_memgrp_);
  //memcpy(dest, static_cast<const void*>(ints), this->blksize());
  //mem_->release_writeonly(dest, offset, blksize_memgrp_);
  if (classdebug() > 0) {
    std::ostringstream oss;
    oss << "i=" << i << ",j=" << j << ": wrote " << blksize_memgrp_ << " bytes to proc=" << ij_proc(i,j) << " offset=" << offset << " buf[0]=" << ints[0] << std::endl;
    ExEnv::outn() << oss.str();
    ExEnv::outn() << scprintf("offsets = %d %d\n",mem_->offset(0), mem_->offset(1));
  }
}

void
DistArray4_MemoryGrp::store_pair_subblock(int i, int j, tbint_type oper_type,
                                          int xstart, int xfence, int ystart, int yfence,
                                          const double *buf)
{
  MPQC_ASSERT(this->active());  //make sure we are active
  MPQC_ASSERT(i>=0 && i<ni());
  MPQC_ASSERT(j>=0 && j<nj());
  MPQC_ASSERT(is_local(i,j));   // store blocks local to this node ONLY

  const bool contiguous = (ystart == 0) && (yfence == ny());
  const int xsize = xfence - xstart;
  const int ysize = yfence - ystart;
  const int xysize = xsize * ysize;
  const int bufsize = xysize * sizeof(double);

  const int ij = ij_index(i,j);
  const int local_ij_index = ij / ntasks();
  const double* srcbuf = buf;
  double* destbuf = static_cast<double*>(mem_->localdata()) + blksize_memgrp_*(num_te_types()*local_ij_index + oper_type)/sizeof(double) + (xstart * ny() + ystart);

  // copy integrals to MemoryGrp
  if (contiguous) { // write all at once
    std::copy(srcbuf, srcbuf + xysize, destbuf);
  }
  else {
    for(int x=0; x<xsize; ++x, destbuf+=ny(), srcbuf+=ysize) {
      std::copy(srcbuf, srcbuf + ysize, destbuf);
    }
  }
}

void
DistArray4_MemoryGrp::deactivate()
{
  DistArray4::deactivate();
  mem_->sync();
}

const double *
DistArray4_MemoryGrp::retrieve_pair_block(int i, int j, tbint_type oper_type,
                                          double* buf) const
{
  MPQC_ASSERT(this->active());  //make sure we are active
  MPQC_ASSERT(i>=0 && i<ni());
  MPQC_ASSERT(j>=0 && j<nj());
  const int ij = ij_index(i,j);
  struct PairBlkInfo *pb = &pairblk_[ij];
  // if it's local all is easy
  // if it's not -- read it if it's not already available;
  if (!is_local(i,j) && pb->ints_[oper_type] == 0) {
    const int local_ij_index = ij / ntasks();
    // this offset could not be computed in ::init(), compute here and store so that release can reuse its value
    pairblk_[ij].offset_ = mem_->offset(ij_proc(i,j)) + (distsize_t)local_ij_index*blksize_memgrp_*num_te_types();
    const distsize_t offset = pb->offset_ + (distsize_t)oper_type*blksize_memgrp_;
    if (classdebug() > 0) {
      std::ostringstream oss;
      oss << indent << "retrieving remote block:  i,j=" << i << "," << j << " oper_type=" << oper_type << " offset=" << offset << endl;
      ExEnv::outn() << oss.str();
    }
    pb->ints_[oper_type] = (double *) mem_->obtain_readonly(offset, blksize());
  }
  pb->refcount_[oper_type] += 1;
  if (classdebug() > 0) {
    std::ostringstream oss;
    oss << indent << me() << ":refcount=" << pb->refcount_[oper_type]
                  << ": i = " << i << " j = " << j << " tbint_type = " << oper_type  << " ptr = " << pb->ints_[oper_type] << endl;
    ExEnv::outn() << oss.str();
  }

  if (buf) {
    // if user provided the buffer, copy to it
    std::copy(pb->ints_[oper_type], pb->ints_[oper_type] + this->nxy(), buf);
    return const_cast<const double*>(buf);
  }
  else
    return pb->ints_[oper_type];
}

void
DistArray4_MemoryGrp::retrieve_pair_subblock(int i, int j, tbint_type oper_type,
                                             int xstart, int xfence, int ystart, int yfence,
                                             double* buf) const
{
  MPQC_ASSERT(this->active());  //make sure we are active
  MPQC_ASSERT(i>=0 && i<ni());
  MPQC_ASSERT(j>=0 && j<nj());
  const bool contiguous = (ystart == 0) && (yfence == ny());
  const int xsize = xfence - xstart;
  const int ysize = yfence - ystart;
  const int xysize = xsize * ysize;
  const int bufsize = xysize * sizeof(double);

  const int ij = ij_index(i,j);
  struct PairBlkInfo *pb = &pairblk_[ij];
  // if it's not local nor available, read
  if (!is_local(i,j) && pb->ints_[oper_type] == 0) {
    const int local_ij_index = ij / ntasks();
    const distsize_t ijxy_offset = mem_->offset(ij_proc(i,j))
                                  + (distsize_t)local_ij_index*blksize_memgrp_*num_te_types()
                                  + (distsize_t)oper_type*blksize_memgrp_
                                  + (distsize_t)(xstart*ny() + ystart)*sizeof(double);

    if (contiguous) { // read all at once
      double* tmpbuf = (double *) mem_->obtain_readonly(ijxy_offset, bufsize);
      std::copy(tmpbuf, tmpbuf + xysize, buf);
      mem_->release_readonly(tmpbuf, ijxy_offset, bufsize);
    }
    else { // read row by row
      distsize_t offset = ijxy_offset;
      double* outbuf = buf;
      const distsize_t stride = ny()*sizeof(double);
      const int rdsize = ysize*sizeof(double);
      for(int x=0; x<xsize; ++x, offset+=stride, outbuf+=ysize) {
        double* tmpbuf = (double *) mem_->obtain_readonly(offset, rdsize);
        std::copy(tmpbuf, tmpbuf + ysize, outbuf);
        mem_->release_readonly(tmpbuf, offset, rdsize);
      }
    }
  }
  else { // it's local or available
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
DistArray4_MemoryGrp::release_pair_block(int i, int j, tbint_type oper_type) const
{
  MPQC_ASSERT(this->active());  //make sure we are active
  int ij = ij_index(i,j);
  struct PairBlkInfo *pb = &pairblk_[ij];
  if (pb->refcount_[oper_type] <= 0) {
    ExEnv::outn() << indent << me() << ":refcount=0: i = " << i << " j = " << j << " tbint_type = " << oper_type << endl;
    throw ProgrammingError("Logic error: DistArray4_MemoryGrp::release_pair_block: refcount is already zero!",__FILE__,__LINE__);
  }
  if (classdebug() > 0)
    ExEnv::outn() << indent << me() << ":refcount=" << pb->refcount_[oper_type]
                  << ": i = " << i << " j = " << j << " tbint_type = " << oper_type << " ptr = " << pb->ints_[oper_type] << " ptr[0] = " << pb->ints_[oper_type][0] << endl;
  if (!is_local(i,j) && pb->ints_[oper_type] != NULL && pb->refcount_[oper_type] == 1) {
    mem_->release_readonly(const_cast<void*>(reinterpret_cast<const void*>(pb->ints_[oper_type])),pb->offset_ + oper_type*blksize_memgrp_,blksize());
    pb->offset_ = ULONG_MAX;
    pb->ints_[oper_type] = NULL;
  }
  pb->refcount_[oper_type] -= 1;
}

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
