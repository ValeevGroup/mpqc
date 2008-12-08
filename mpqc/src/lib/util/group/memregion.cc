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

#ifdef __GNUC__
#pragma implementation
#endif

#include <cassert>
#include <util/group/memregion.h>

using namespace sc;

static ClassDesc MemoryGrpRegion_cd(
  typeid(MemoryGrpRegion),"MemoryGrpRegion",1,"public MemoryGrp",
  0, 0, 0);

MemoryGrpRegion::MemoryGrpRegion(const Ref<MemoryGrp>& host, size_t mem_offset, size_t max_size) :
  MemoryGrp(), reserve_(mem_offset,max_size)
{
  n_ = host->n();
  me_ = host->me();
}

MemoryGrpRegion::~MemoryGrpRegion()
{
}

void
MemoryGrpRegion::set_localsize(size_t localsize)
{
  delete[] offsets_;
  offsets_ = new distsize_t[n_ + 1];
  assert(false);
}

void*
MemoryGrpRegion::localdata()
{
  return static_cast<void*>( static_cast<char*>(host_->localdata()) + reserve_.start()) ;
}

void *
MemoryGrpRegion::obtain_readwrite(distsize_t offset, int size)
{
  obtain_local_lock(offset-localoffset(), offset-localoffset()+size);
//  return &data_[distsize_to_size(offset)];
}

void *
MemoryGrpRegion::obtain_readonly(distsize_t offset, int size)
{
//  return &data_[distsize_to_size(offset)];
}

void *
MemoryGrpRegion::obtain_writeonly(distsize_t offset, int size)
{
//  return &data_[distsize_to_size(offset)];
}

void
MemoryGrpRegion::release_readonly(void *data, distsize_t offset, int size)
{
}

void
MemoryGrpRegion::release_writeonly(void *data, distsize_t offset, int size)
{
}

void
MemoryGrpRegion::release_readwrite(void *data, distsize_t offset, int size)
{
  release_local_lock(offset-localoffset(), offset-localoffset()+size);
}

void
MemoryGrpRegion::sync()
{
}

#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
