//
// memproc.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
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

#ifndef _util_group_memproc_cc
#define _util_group_memproc_cc

#ifdef __GNUC__
#pragma implementation
#endif

#include <util/group/memproc.h>

using namespace sc;

static ClassDesc ProcMemoryGrp_cd(
  typeid(ProcMemoryGrp),"ProcMemoryGrp",1,"public MemoryGrp",
  create<ProcMemoryGrp>, create<ProcMemoryGrp>, 0);

ProcMemoryGrp::ProcMemoryGrp()
{
  data_ = 0;
  offsets_ = 0;
}

ProcMemoryGrp::ProcMemoryGrp(const Ref<KeyVal>& keyval):
  MemoryGrp(keyval)
{
  data_ = 0;
  offsets_ = 0;
}

ProcMemoryGrp::~ProcMemoryGrp()
{
  delete[] data_;
}

void
ProcMemoryGrp::set_localsize(size_t localsize)
{
  delete[] offsets_;
  delete[] data_;
  offsets_ = new distsize_t[2];
  offsets_[0] = 0;
  offsets_[1] = localsize;
  n_ = 1;
  me_ = 0;
  data_ = new char[localsize];
}

void *
ProcMemoryGrp::localdata()
{
  return data_;
}

void *
ProcMemoryGrp::obtain_readwrite(distsize_t offset, int size)
{
  obtain_local_lock(offset-localoffset(), offset-localoffset()+size);
  return &data_[distsize_to_size(offset)];
}

void *
ProcMemoryGrp::obtain_readonly(distsize_t offset, int size)
{
  return &data_[distsize_to_size(offset)];
}

void *
ProcMemoryGrp::obtain_writeonly(distsize_t offset, int size)
{
  return &data_[distsize_to_size(offset)];
}

void
ProcMemoryGrp::release_readonly(void *data, distsize_t offset, int size)
{
}

void
ProcMemoryGrp::release_writeonly(void *data, distsize_t offset, int size)
{
}

void
ProcMemoryGrp::release_readwrite(void *data, distsize_t offset, int size)
{
  release_local_lock(offset-localoffset(), offset-localoffset()+size);
}

void
ProcMemoryGrp::sync()
{
}

#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
