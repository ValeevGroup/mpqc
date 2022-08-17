//
// memregion.cc
//
// Copyright (C) 2008 Edward Valeev
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

#ifndef _util_group_memregion_cc
#define _util_group_memregion_cc

#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <util/misc/scexception.h>
#include <util/group/memregion.h>
#include <util/group/message.h>

using namespace sc;

#ifdef HAVE_HRECV
#  define DISABLE do { masktrap(1); ExEnv::outn().flush(); } while(0)
#  define ENABLE do { ExEnv::outn().flush(); masktrap(0); } while(0)
   extern "C" {
       long masktrap(long state);
     }
#else
#  define DISABLE ExEnv::outn().flush()
#  define ENABLE ExEnv::outn().flush()
#endif

#define PRINTF(args) do { DISABLE; \
                          ExEnv::outn() << scprintf args ; \
                          ExEnv::outn().flush(); \
                          ENABLE; \
                         } while(0)

// comment out these two lines to produce diagnostic output
#undef PRINTF
#define PRINTF(args)

static ClassDesc MemoryGrpRegion_cd(
  typeid(MemoryGrpRegion),"MemoryGrpRegion",1,"public MemoryGrp",
  0, 0, 0);

MemoryGrpRegion::MemoryGrpRegion(const Ref<MemoryGrp>& host, size_t host_offset, size_t max_size) :
  MemoryGrp(), host_(host), reserve_(host_offset,max_size), host_offsets_(host_->n(),0)
{
  if (only_allow_double_aligned_regions_ &&
      (reserve_.start()%sizeof(double) != 0 ||
       reserve_.size()%sizeof(double) != 0)
     )
    throw ProgrammingError("MemoryGrpRegion::MemoryGrpRegion -- only sizeof(double)-aligned regions are supported",__FILE__,__LINE__);

  // host must be initialized already and have enough memory
  if (reserve_.start() + reserve_.size() > host->localsize()) {
    std::ostringstream oss;
    oss << "MemoryGrpRegion: not enough memory, got free region at " << reserve_.start()
        << ", need " << reserve_.size() << " bytes, localsize = " << host->localsize() << std::endl;
    throw std::runtime_error(oss.str().c_str());
  }
  // add region to the map
  map_.insert(host_,reserve_);

  init();
}

MemoryGrpRegion::MemoryGrpRegion(const Ref<MemoryGrp>& host, size_t max_size) :
  MemoryGrp(), host_(host),
  // find_free finds a free segment of size max_size.
  reserve_(map_.find_free(host_,max_size),max_size),
  host_offsets_(host_->n(),0)
{
  if (only_allow_double_aligned_regions_ &&
      (reserve_.start()%sizeof(double) != 0 ||
       reserve_.size()%sizeof(double) != 0)
     )
    throw ProgrammingError("MemoryGrpRegion::MemoryGrpRegion -- only sizeof(double)-aligned regions are supported",__FILE__,__LINE__);

  // host must be initialized already and have enough memory
  if (reserve_.start() + reserve_.size() > host->localsize()) {
    std::ostringstream oss;
    oss << "MemoryGrpRegion: not enough memory, got free region at " << reserve_.start()
        << ", need " << reserve_.size() << " bytes, localsize = " << host->localsize() << std::endl;
    throw std::runtime_error(oss.str().c_str());
  }
  // add region to the map
  map_.insert(host_,reserve_);

  init();
}

void
MemoryGrpRegion::init()
{
  // replicate host_offset
  const int nnodes = host_->n();
  long int* host_offsets = new long int[nnodes];
  for(int node=0; node<nnodes; ++node)
    host_offsets[node] = 0;
  host_offsets[host_->me()] = reserve_.start();
  Ref<MessageGrp> msg = MessageGrp::get_default_messagegrp();
  msg->sum(host_offsets,nnodes);

  // compute host_offsets_
  host_offsets_.resize(nnodes);
  for(int node=0; node<nnodes; ++node)
    host_offsets_[node] = host_->offset(node) + host_offsets[node];

  // MemoryGrp constructor doesn't do much, initialize its protected members here
  n_ = host_->n();
  me_ = host_->me();
  // initialize MemoryGrp::offsets_ -- it has nnodes+1 elements!
  offsets_ = new distsize_t[nnodes+1];
  for(int node=0; node<=nnodes; ++node) offsets_[node] = 0;

  // cleanup
  delete[] host_offsets;

  // print debug info
  if (classdebug()) {
    this->print(ExEnv::out0());
  }

}

MemoryGrpRegion::~MemoryGrpRegion()
{
  map_.erase(host_,reserve_);
}

void
MemoryGrpRegion::set_localsize(size_t localsize)
{
  MPQC_ASSERT(localsize <= host_->localsize());
  MPQC_ASSERT(localsize <= reserve_.size());

  // replicate sizes
  const int nnodes = host_->n();
  long int* sizes = new long int[nnodes];
  for(int node=0; node<nnodes; ++node) sizes[node] = 0;
  sizes[host_->me()] = localsize;
  Ref<MessageGrp> msg = MessageGrp::get_default_messagegrp();
  msg->sum(sizes,nnodes);

  // compute offsets_
  offsets_[0] = 0;
  // remember, offsets_ is nnodes+1 long
  for(int node=1; node<=nnodes; ++node)
    offsets_[node] = offsets_[node-1] + sizes[node-1];

  // cleanup
  delete[] sizes;

  if (classdebug()) {
    std::ostringstream oss;
    oss << "MemoryGrpRegion::set_localsize(" << localsize << ") called" << std::endl;
    ExEnv::out0() << indent << oss.str() << incindent;
    this->print(ExEnv::out0());
    ExEnv::out0() << decindent;
  }
}

void*
MemoryGrpRegion::localdata()
{
  return static_cast<void*>( static_cast<char*>(host_->localdata()) + reserve_.start()) ;
}

distsize_t
MemoryGrpRegion::offset_to_host_offset(const distsize_t& offset) const
{
  // find the last node with host_offset less than or equal to the requested offset
  int node = host_->n() - 1;
  for(; node>=0; --node)
    if (offsets_[node] <= offset)
      break;
  const distsize_t host_offset = host_offsets_[node] + (offset - offsets_[node]);
  if (classdebug() > 1) {
    PRINTF(("MemoryGrpRegion::offset_to_host_offset -- offset = %ld host_offset = %ld\n", (size_t)offset, (size_t)host_offset));
  }
  return host_offset;
}

void *
MemoryGrpRegion::obtain_readwrite(distsize_t offset, size_t size)
{
  const distsize_t host_offset = offset_to_host_offset(offset);
  return host_->obtain_readwrite(host_offset,size);
}

void *
MemoryGrpRegion::obtain_readonly(distsize_t offset, size_t size)
{
  const distsize_t host_offset = offset_to_host_offset(offset);
  return host_->obtain_readonly(host_offset,size);
}

void *
MemoryGrpRegion::obtain_writeonly(distsize_t offset, size_t size)
{
  const distsize_t host_offset = offset_to_host_offset(offset);
  return host_->obtain_writeonly(host_offset,size);
}

void
MemoryGrpRegion::release_readonly(void *data, distsize_t offset, size_t size)
{
  const distsize_t host_offset = offset_to_host_offset(offset);
  host_->release_readonly(data,host_offset,size);
}

void
MemoryGrpRegion::release_writeonly(void *data, distsize_t offset, size_t size)
{
  const distsize_t host_offset = offset_to_host_offset(offset);
  return host_->release_writeonly(data,host_offset,size);
}

void
MemoryGrpRegion::release_readwrite(void *data, distsize_t offset, size_t size)
{
  const distsize_t host_offset = offset_to_host_offset(offset);
  return host_->release_readwrite(data,host_offset,size);
}

void
MemoryGrpRegion::sum_reduction_on_node(double *data, size_t doffset, size_t dsize, int node)
{
  const distsize_t offset = doffset * sizeof(double);
  const distsize_t host_offset = offset_to_host_offset(offset);
  host_->sum_reduction_on_node(data,host_offset/sizeof(double),dsize,node);
}

#if 0
void
MemoryGrpRegion::sum_reduction(double *data, distsize_t doffset, size_t dsize)
{
  // doffset is in terms of doubles
  const distsize_t offset = doffset * sizeof(double);
  const size_t size = dsize * sizeof(double);

  if (offset + size > totalsize()) {
    std::ostrinstream oss;
    oss << "MemoryGrpRegion::sum_reduction: arg out of range:"
      << " offset = " << offset
      << " length = " << length
      << " totalsize() = " << totalsize()
      << std::endl;
    throw ProgrammingError(oss.str().c_str(),__FILE__,__LINE__);
  }

  double *source_data = (double*) obtain_readwrite(offset, size);
  for (size_t i=0; i<dsize; i++) {
    source_data[i] += data[i];
  }
  release_readwrite((void*) source_data, offset, size);
}
#endif

void
MemoryGrpRegion::sync()
{
  host_->sync();
}

void
MemoryGrpRegion::catchup()
{
  host_->catchup();
}

Ref<MemoryGrp>
MemoryGrpRegion::clone()
{
  Ref<MemoryGrp> result = new MemoryGrpRegion(host_, reserve_.size());
  return result;
}

void
MemoryGrpRegion::activate()
{
  host_->activate();
}

void
MemoryGrpRegion::deactivate()
{
  host_->deactivate();
}

void *
MemoryGrpRegion::malloc_local(size_t nbytes) {
  return host_->malloc_local(nbytes);
}

void
MemoryGrpRegion::free_local(void*& data) {
  host_->free_local(data);
}

void
MemoryGrpRegion::print(std::ostream& o) const
{
  std::ostringstream oss;
  oss << scprintf("MemoryRegionGrp (node %d):\n", me());
  oss << scprintf("%d: n = %d max_size = %ld\n", me(), n(), reserve_.size());
  for (int i=0; i<=n_; i++) {
    oss << me() << ": offset[" << i << "] = " << std::setw(10) << offsets_[i];
    if (i < n_)
      oss << "  host_offset[" << i << "] = " << std::setw(10) << host_offsets_[i];
    oss << std::endl;
  }
  oss << "Host MemoryGrp (class " << host_->class_name() << ")" << std::endl << incindent;
  host_->print(oss);
  o << oss.str() << decindent;
}

/////////////////////////////////////////////////////////////////////////////

MemoryGrpRegion::HostToRegionsMap MemoryGrpRegion::map_;

MemoryGrpRegion::HostToRegionsMap::HostToRegionsMap() :
 lock_(ThreadGrp::get_default_threadgrp()->new_lock())
{
}

MemoryGrpRegion::HostToRegionsMap::~HostToRegionsMap() {
  typedef ImplType::iterator iter;
  for(iter h=impl_.begin();
      h!=impl_.end();
      ++h) {
    delete h->second;
  }
}

// find the regions list for this host
const MemoryGrpRegion::LocalRegions*
MemoryGrpRegion::HostToRegionsMap::regions(const MemoryGrp* host) const {
  // assumed already locked
  const LocalRegions* result = 0;
  ImplType::const_iterator regions_iter = impl_.find(host);
  if (regions_iter != impl_.end())
    result = regions_iter->second;
  return result;
}

void
MemoryGrpRegion::HostToRegionsMap::insert(const MemoryGrp* host, const LocalRegion& region)
{

  // one thread at a time
  lock_->lock();

  LocalRegions* regions_ptr = const_cast<LocalRegions*>(this->regions(host));
  if (regions_ptr) {
    LocalRegions& regions = *regions_ptr;
    typedef LocalRegions::iterator iter;
    typedef LocalRegions::reverse_iterator riter;
    // find region V which should follow region I
    iter v = std::find_if(regions.begin(),regions.end(),std::bind(std::greater<LocalRegion>(),std::placeholders::_1,region));
    if (v != regions.end()) {
      // if found -- check if the two regions overlap ...
      if (region.start() + region.size() > v->start()) {
        lock_->unlock();
        throw ProgrammingError("MemoryGrpRegion::HostToRegionsMap::insert -- region overlaps with its neighbor",__FILE__,__LINE__);
      }
      // and check if the region to precede this region exists
      // if yes, make sure they don't overlap and insert
      if (v != regions.begin()) {
        iter vprev = v; --vprev;
        if (vprev->start() + vprev->size() > region.start()) {
          lock_->unlock();
          throw ProgrammingError("MemoryGrpRegion::HostToRegionsMap::insert -- region overlaps with its neighbor",__FILE__,__LINE__);
        }
        regions.insert(v,region);
      }
      else { // or simply prepend
        regions.push_front(region);
      }
    } // v != regions.rend()
    else { // else check if the region to precede this region doesn't overlap, then append
      riter vprev = regions.rbegin();
      if (vprev->start() + vprev->size() > region.start()) {
        lock_->unlock();
        throw ProgrammingError("MemoryGrpRegion::HostToRegionsMap::insert -- region overlaps with its neighbor",__FILE__,__LINE__);
      }
      regions.push_back(region);
    }
  }
  else {
    LocalRegions* regions_ptr = new LocalRegions;
    regions_ptr->push_back(region);
    impl_[host] = regions_ptr;
  }

  lock_->unlock();
}

void
MemoryGrpRegion::HostToRegionsMap::erase(const MemoryGrp* host, const LocalRegion& region) {

  // one thread at a time
  lock_->lock();

  LocalRegions* regions_ptr = const_cast<LocalRegions*>(this->regions(host));
  if (regions_ptr == 0) {
    lock_->unlock();
    throw ProgrammingError("MemoryGrpRegion::HostToRegionsMap::erase -- no regions associated with this host",__FILE__,__LINE__);
  }
  LocalRegions& regions = *regions_ptr;

  typedef LocalRegions::iterator iter;
  // find region V which should follow region I
  iter v = std::find(regions.begin(),regions.end(),region);
  if (v != regions.end())
    regions.erase(v);
  else { // else not found -- throw
    lock_->unlock();
    throw ProgrammingError("MemoryGrpRegion::HostToRegionsMap::erase -- this region not found",__FILE__,__LINE__);
  }

  // remove regions, if empty
  if (regions.empty())
    impl_.erase(host);

  lock_->unlock();
}

size_t
MemoryGrpRegion::HostToRegionsMap::find_free(const MemoryGrp* host, size_t size) const
{
  // one thread at a time
    lock_->lock();

    size_t result = 0;
    LocalRegions* regions_ptr = const_cast<LocalRegions*>(this->regions(host));
    // if regions_ptr == 0, return 0
    if (regions_ptr != 0) {
      LocalRegions& regions = *regions_ptr;
      // regions are ordered, find the first gap that is big enough, otherwise assume free space past the last block
      typedef LocalRegions::iterator iter;
      const iter end = regions.end();
      for(iter block=regions.begin(); block!=end; ++block) {
        const size_t gap_start = block->start() + block->size();
        result = gap_start;
        // if the gap between the end of this block and the next one if large enough, stop
        iter next_block = block; ++next_block;
        if (next_block != end) {
          const size_t gap_size = next_block->start() - gap_start;
          if (gap_size >= size) {
            break;
          }
        }
      }
    }

    lock_->unlock();
    return result;
}

#endif // header guard

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
