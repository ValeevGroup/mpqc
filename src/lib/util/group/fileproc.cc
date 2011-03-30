//
// fileproc.cc
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

#ifndef _util_group_fileproc_cc
#define _util_group_fileproc_cc

#include <util/group/fileproc.h>

using namespace sc;

static ClassDesc ProcFileGrp_cd(
  typeid(ProcFileGrp),"ProcFileGrp",1,"public FileGrp",
  0, create<ProcFileGrp>, 0);

ProcFileGrp::ProcFileGrp()
{
  data_ = 0;
  offsets_ = 0;
}

ProcFileGrp::ProcFileGrp(const Ref<KeyVal>& keyval):
  FileGrp(keyval)
{
  data_ = 0;
  offsets_ = 0;
}

ProcFileGrp::~ProcFileGrp()
{
  delete[] data_;
}

ProcFileGrp*
ProcFileGrp::clone()
{
  return new ProcFileGrp;
}

void
ProcFileGrp::set_localsize(size_t localsize)
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
ProcFileGrp::localdata()
{
  return data_;
}

void *
ProcFileGrp::obtain_readwrite(distsize_t offset, int size)
{
  obtain_local_lock(offset-localoffset(), offset-localoffset()+size);
  return &data_[distsize_to_size(offset)];
}

void *
ProcFileGrp::obtain_readonly(distsize_t offset, int size)
{
  return &data_[distsize_to_size(offset)];
}

void *
ProcFileGrp::obtain_writeonly(distsize_t offset, int size)
{
  return &data_[distsize_to_size(offset)];
}

void
ProcFileGrp::release_readonly(void *data, distsize_t offset, int size)
{
}

void
ProcFileGrp::release_writeonly(void *data, distsize_t offset, int size)
{
}

void
ProcFileGrp::release_readwrite(void *data, distsize_t offset, int size)
{
  release_local_lock(offset-localoffset(), offset-localoffset()+size);
}

void
ProcFileGrp::sync()
{
}

#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
