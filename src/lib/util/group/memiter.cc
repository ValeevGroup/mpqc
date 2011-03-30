//
// memiter.cc
//
// derived from memamsg.cc
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

#ifndef _util_group_memiter_cc
#define _util_group_memiter_cc

#include <util/group/memiter.h>

using namespace sc;

///////////////////////////////////////////////////////////////////////
// The MemoryIter class

int
MemoryIter::local(distsize_t offset, int size, int node)
{
  if (offset >= offsets_[node] && offset + size <= offsets_[node+1])
      return 1;
  return 0;
}

MemoryIter::MemoryIter(void *data,
                       distsize_t *offsets,
                       int n):
  offsets_(offsets),
  n_(n),
  data_(data),
  ready_(0)
{
}

void
MemoryIter::begin(distsize_t offset, int size)
{
  offset_ = offset;
  size_ = size;

  current_data_ = (char *) data_;

  distsize_t fence = offset + size;

  for (node_ = 0; node_ < n_; node_++) {
      if (offset_ < offsets_[node_ + 1]) {
          current_offset_ = distsize_to_size(offset_ - offsets_[node_]);
          if (fence <= offsets_[node_ + 1]) {
              current_size_ = size_;
            }
          else {
              current_size_ = distsize_to_size(offsets_[node_ + 1] - offset_);
            }
          ready_ = 1;
          return;
        }
    }

  // couldn't find the requested data, this is probably from an
  // invalid argument
  ready_ = 0;
}

void
MemoryIter::next()
{
  distsize_t fence = offset_ + size_;

  if (fence <= offsets_[node_ + 1]) {
      ready_ = 0;
    }
  else {
      node_++;
      current_data_ += current_size_;
      if (fence <= offsets_[node_ + 1]) {
          current_size_ = size_ - distsize_to_size(offsets_[node_] - offset_);
        }
      else {
          current_size_ = distsize_to_size(offsets_[node_+1]-offsets_[node_]);
        }
      current_offset_ = 0;
    }
}

#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
