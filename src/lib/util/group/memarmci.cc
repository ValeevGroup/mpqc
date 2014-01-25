//
// memarmci.cc
// based on memshm.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
// Maintainer: SNL
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

#ifndef _util_group_memarmci_cc
#define _util_group_memarmci_cc

// ARMCI may include MPI/C++ headers
#define OMPI_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
extern "C" {
#include <armci.h>
}

#include <stdexcept>

#include <util/misc/formio.h>
#include <util/misc/consumableresources.h>
#include <util/misc/scexception.h>
#include <util/group/memarmci.h>

using namespace sc;

static ClassDesc ARMCIMemoryGrp_cd(
  typeid(ARMCIMemoryGrp),"ARMCIMemoryGrp",1,"public RDMAMemoryGrp",
  0, create<ARMCIMemoryGrp>, 0);

ARMCIMemoryGrp::ARMCIMemoryGrp(const Ref<MessageGrp>& msg):
  RDMAMemoryGrp(msg)
{
  init();
}

ARMCIMemoryGrp::ARMCIMemoryGrp(const Ref<KeyVal>& keyval):
  RDMAMemoryGrp(keyval)
{
  init();
}

void
ARMCIMemoryGrp::init()
{
  armci_lock_ = ThreadGrp::get_default_threadgrp()->new_lock();
  //debug_ = 1;
  all_data_ = 0;
  ARMCI_Init();
}

void
ARMCIMemoryGrp::finalize()
{
  set_localsize(0);
  ARMCI_Finalize();
}

void
ARMCIMemoryGrp::set_localsize(size_t localsize)
{
  ARMCI_AllFence();

  const size_t current_localsize = offsets_[me()+1] - offsets_[me()];

  // this will initialize the offsets_ array
  RDMAMemoryGrp::set_localsize(localsize);

  if (all_data_) {
      ARMCI_Free(data_);
      unmanage_array(data_);
      delete[] all_data_;
      all_data_ = 0;
      data_ = 0;
      ARMCI_Destroy_mutexes();
    }

  if (localsize == 0) return;

  all_data_ = new void*[n()];
  int r;
  r = ARMCI_Malloc(all_data_, localsize);
  data_ = reinterpret_cast<char*>(all_data_[me()]);
  manage_array(data_, localsize);

  if (debug_) {
    for (int i=0; i<n(); i++) {
      std::cout << me() << ": all_data[" << i
		<< "] = " << all_data_[i] << std::endl;
    }
  }

  ARMCI_Create_mutexes(1);
}

void
ARMCIMemoryGrp::retrieve_data(void *data, int node, long offset,
                              long size, int lock)
{
  if (armci_lock_) armci_lock_->lock();
  if (lock) ARMCI_Lock(0, node);
  ARMCI_Get(reinterpret_cast<char*>(all_data_[node])+offset, data, size, node);
  if (armci_lock_) armci_lock_->unlock();
}

void
ARMCIMemoryGrp::replace_data(void *data, int node, long offset,
                             long size, int unlock)
{
  if (armci_lock_) armci_lock_->lock();
  ARMCI_Put(data, reinterpret_cast<char*>(all_data_[node])+offset, size, node);
  if (unlock) {
      ARMCI_Fence(node);
      ARMCI_Unlock(0, node);
    }
  if (armci_lock_) armci_lock_->unlock();
}

void
ARMCIMemoryGrp::sum_data(double *data, int node, long offset, long size)
{
  long doffset = offset/sizeof(double);
  long dsize = size/sizeof(double);

  void *src = data;
  void *dst = reinterpret_cast<double*>(all_data_[node])+doffset;

  armci_giov_t acc_dat;
  acc_dat.src_ptr_array = &src;
  acc_dat.dst_ptr_array = &dst;
  acc_dat.bytes = dsize * sizeof(double);
  acc_dat.ptr_array_len = 1;
  double scale = 1.0;

  if (debug_) {
      std::cout << me() << ": summing " << dsize
                << " doubles from "
                << (void*)src
                << " to "
                << (void*)dst
                << " on " << node
                << " (base dest=" << (void*)all_data_[node] << ")"
                << std::endl;
      for (int i=0; i<dsize; i++) {
          std::cout << me() << ": src[" << i << "] = "
                    << data[i] << std::endl;
        }
//        for (int i=0; i<dsize; i++) {
//            std::cout << me() << ": dst[" << i << "] = "
//                      << ((double*)(all_data_[node]))[doffset+i]
//                      << std::endl;
//          }
    }

  if (armci_lock_) armci_lock_->lock();
  // Original code sending all data at once:
  // ARMCI_AccV(ARMCI_ACC_DBL, &scale, &acc_dat, 1, node);
  // Hack to send smaller chunks to not overflow buffers in ARMCI:
  int incr = 32768;
  for (int i=0; i<size; i+=incr) {
      void *tsrc = (&(((char*)src)[i]));
      void *tdst = (&(((char*)dst)[i]));
      acc_dat.src_ptr_array = &tsrc;
      acc_dat.dst_ptr_array = &tdst;
      if (size - i > incr) acc_dat.bytes = incr;
      else acc_dat.bytes = (size-i);
      acc_dat.ptr_array_len = 1;
      ARMCI_AccV(ARMCI_ACC_DBL, &scale, &acc_dat, 1, node);
    }
  // Send data all at once using the contiguous routine (which does not exist):
  // ARMCI_Acc(ARMCI_ACC_DBL, &scale, src, dst, size, node);
  if (armci_lock_) armci_lock_->unlock();
}

void
ARMCIMemoryGrp::sync()
{
  ARMCI_Barrier();
}

void
ARMCIMemoryGrp::deactivate()
{
  // Really, this is still active after deactivate is called.
  // However, we'll at least make sure that all outstanding
  // requests are finished.
  ARMCI_AllFence();
}

void*
ARMCIMemoryGrp::malloc_local(size_t nbyte)
{
  void* buf = ARMCI_Malloc_local(nbyte);
  if (buf == NULL)
    throw MemAllocFailed("malloc_local -- failed to allocate memory",
                         __FILE__, __LINE__, nbyte, this->class_desc());
  manage_array(buf, nbyte);
  return buf;
}

void
ARMCIMemoryGrp::free_local(void* &data)
{
  ARMCI_Free_local(data);
  unmanage_array(data);
}

ARMCIMemoryGrp::~ARMCIMemoryGrp()
{
  finalize();
}

void
ARMCIMemoryGrp::print(std::ostream &o) const
{
  RDMAMemoryGrp::print(o);
}

#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
