//
// memproc.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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

#define CLASSNAME ProcMemoryGrp
#define HAVE_KEYVAL_CTOR
#define PARENTS public MemoryGrp
#include <util/class/classi.h>
void *
ProcMemoryGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  MemoryGrp::_castdown(cd);
  return do_castdowns(casts,cd);
}

ProcMemoryGrp::ProcMemoryGrp()
{
  data_ = 0;
  offsets_ = 0;
}

ProcMemoryGrp::ProcMemoryGrp(const RefKeyVal& keyval):
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
ProcMemoryGrp::set_localsize(int localsize)
{
  delete[] offsets_;
  delete[] data_;
  offsets_ = new int[2];
  offsets_[0] = 0;
  offsets_[1] = localsize;
  n_ = 1;
  me_ = 0;
  data_ = new char[localsize];
}

void *
ProcMemoryGrp::obtain_readwrite(int offset, int size)
{
  return &data_[offset];
}

void *
ProcMemoryGrp::obtain_readonly(int offset, int size)
{
  return &data_[offset];
}

void
ProcMemoryGrp::release_read(void *data, int offset, int size)
{
}

void
ProcMemoryGrp::release_write(void *data, int offset, int size)
{
}

void
ProcMemoryGrp::sync()
{
}

#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
