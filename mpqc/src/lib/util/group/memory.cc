//
// memory.cc
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

#ifndef _util_group_memory_cc
#define _util_group_memory_cc

#ifdef __GNUC__
#pragma implementation
#endif

#include <util/group/memory.h>

//////////////////////////////////////////////////////////////////////
// MemoryGrpBuf members

template <class data_t>
MemoryGrpBuf<data_t>::MemoryGrpBuf(const RefMemoryGrp & grp)
{
  grp_ = grp;
  locktype_ = None;
}

template <class data_t>
data_t *
MemoryGrpBuf<data_t>::writeonly(distsize_t offset, int length)
{
  if (locktype_ != None) release();
  data_ = (data_t *) grp_->obtain_writeonly(sizeof(data_t)*offset,
                                            sizeof(data_t)*length);
  offset_ = offset;
  length_ = length;
  locktype_ = Write;
  return data_;
}

template <class data_t>
data_t *
MemoryGrpBuf<data_t>::readwrite(distsize_t offset, int length)
{
  if (locktype_ != None) release();
  data_ = (data_t *) grp_->obtain_readwrite(sizeof(data_t)*offset,
                                            sizeof(data_t)*length);
  offset_ = offset;
  length_ = length;
  locktype_ = Write;
  return data_;
}

template <class data_t>
const data_t *
MemoryGrpBuf<data_t>::readonly(distsize_t offset, int length)
{
  if (locktype_ != None) release();
  data_ = (data_t *) grp_->obtain_readonly(sizeof(data_t)*offset,
                                           sizeof(data_t)*length);
  offset_ = offset;
  length_ = length;
  locktype_ = Read;
  return data_;
}

template <class data_t>
data_t *
MemoryGrpBuf<data_t>::writeonly_on_node(int offset, int length, int node)
{
  if (node == -1) node = grp_->me();
  return writeonly(offset + grp_->offset(node)/sizeof(data_t), length);
}

template <class data_t>
data_t *
MemoryGrpBuf<data_t>::readwrite_on_node(int offset, int length, int node)
{
  if (node == -1) node = grp_->me();
  return readwrite(offset + grp_->offset(node)/sizeof(data_t), length);
}

template <class data_t>
const data_t *
MemoryGrpBuf<data_t>::readonly_on_node(int offset, int length, int node)
{
  if (node == -1) node = grp_->me();
  return readonly(offset + grp_->offset(node)/sizeof(data_t), length);
}

template <class data_t>
void
MemoryGrpBuf<data_t>::release()
{
  if (locktype_ == Write)
      grp_->release_write((data_t *)data_,
                          sizeof(data_t)*offset_, sizeof(data_t)*length_);
  if (locktype_ == Read)
      grp_->release_read(data_, sizeof(data_t)*offset_, sizeof(data_t)*length_);

  locktype_ = None;
}

#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
