//
// memmpi2.cc
//
// derived from code
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@aros.ca.sandia.gov>
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

#ifndef _util_group_memmpi2_cc
#define _util_group_memmpi2_cc

#ifdef __GNUC__
#pragma implementation
#endif

#include <util/group/memmpi2.h>
#include <util/group/memiter.h>
#include <util/misc/formio.h>

using namespace std;

///////////////////////////////////////////////////////////////////////
// The MPI2MemoryGrp class

static ClassDesc MPI2MemoryGrp_cd(
  typeid(MPI2MemoryGrp),"MPI2MemoryGrp",1,"public MsgMemoryGrp",
  0, create<MPI2MemoryGrp>, 0);

MPI2MemoryGrp::MPI2MemoryGrp(const Ref<MessageGrp>& msg):
  MsgMemoryGrp(msg)
{
  rma_win_ = MPI_WIN_NULL;
  data_ = 0;
}

MPI2MemoryGrp::MPI2MemoryGrp(const Ref<KeyVal>& keyval):
  MsgMemoryGrp(keyval)
{
  rma_win_ = MPI_WIN_NULL;
  data_ = 0;
}

MPI2MemoryGrp::~MPI2MemoryGrp()
{
  deactivate();
  delete[] data_;
}

void
MPI2MemoryGrp::set_localsize(size_t localsize)
{
  deactivate();
  MsgMemoryGrp::set_localsize(localsize);
  delete[] data_;
  data_ = new char[localsize];
  activate();
}

void
MPI2MemoryGrp::activate()
{
  MPI_Info info;
  MPI_Info_create(&info);
  MPI_Info_set(info, "no_locks", "true");
  int r = MPI_Win_create(data_, localsize(), sizeof(int), info,
                         MPI_COMM_WORLD, &rma_win_);
  MPI_Info_free(&info);

  if (r != MPI_SUCCESS) {
    ExEnv::err() << "MPI2MemoryGrp::set_localsize("
                 << localsize() << ") failed on "
                 << me() << endl;
    abort();
  }

  sync();
}

void
MPI2MemoryGrp::deactivate()
{
  if (rma_win_ != MPI_WIN_NULL) MPI_Win_free(&rma_win_);
}

void
MPI2MemoryGrp::retrieve_data(void *data, int node, int offset, int size)
{
  int r = MPI_Get(data, size, MPI_BYTE,
                  node, offset, size, MPI_BYTE, rma_win_);
  if (r != MPI_SUCCESS) {
    ExEnv::err() << "MPI2MemoryGrp::retrieve_data failed" << endl;
    abort();
  }
}

void
MPI2MemoryGrp::replace_data(void *data, int node, int offset, int size)
{
  int r = MPI_Put(data, size, MPI_BYTE,
                  node, offset, size, MPI_BYTE, rma_win_);
  if (r != MPI_SUCCESS) {
    ExEnv::err() << "MPI2MemoryGrp::replace_data failed" << endl;
    abort();
  }
}

void
MPI2MemoryGrp::sum_data(double *data, int node, int offset, int size)
{
  int doffset = offset/sizeof(double);
  int dsize = size/sizeof(double);

  int r = MPI_Accumulate(data, dsize, MPI_DOUBLE,
                         node, doffset, dsize, MPI_DOUBLE,
                         MPI_SUM, rma_win_);
  if (r != MPI_SUCCESS) {
    ExEnv::err() << "MPI2MemoryGrp::sum_data failed" << endl;
    abort();
  }
}

void *
MPI2MemoryGrp::obtain_writeonly(distsize_t offset, int size)
{
  void *data = (void *) new char[size];
  return data;
}

void *
MPI2MemoryGrp::obtain_readwrite(distsize_t offset, int size)
{
  static int warn = 1;
  if (warn) {
    warn = 0;
    ExEnv::err() << "MPI2MemoryGrp: readwrite doesn't properly lock"
                 << endl;
  }
  void *data = (void *) new char[size];
  MemoryIter i(data, offsets_, n());
  for (i.begin(offset, size); i.ready(); i.next()) {
      retrieve_data(i.data(), i.node(), i.offset(), i.size());
    }
  return data;
}

void *
MPI2MemoryGrp::obtain_readonly(distsize_t offset, int size)
{
  void *data = (void *) new char[size];
  MemoryIter i(data, offsets_, n());
  for (i.begin(offset, size); i.ready(); i.next()) {
      retrieve_data(i.data(), i.node(), i.offset(), i.size());
    }
  return data;
}

void
MPI2MemoryGrp::sum_reduction(double *data, distsize_t doffset, int dsize)
{
  distsize_t offset = doffset * sizeof(double);
  int size = dsize * sizeof(double);
  MemoryIter i(data, offsets_, n());
  for (i.begin(offset, size); i.ready(); i.next()) {
      sum_data((double*)i.data(), i.node(), i.offset(), i.size());
    }
}

void
MPI2MemoryGrp::sum_reduction_on_node(double *data, size_t doffset,
                                     int dlength, int node)
{
  if (node == -1) node = me();

  sum_data(data, node, sizeof(double)*doffset, sizeof(double)*dlength);
}

void
MPI2MemoryGrp::release_readonly(void *data, distsize_t offset, int size)
{
  delete[] data;
}

void
MPI2MemoryGrp::release_writeonly(void *data, distsize_t offset, int size)
{
  MemoryIter i(data, offsets_, n());
  for (i.begin(offset, size); i.ready(); i.next()) {
      replace_data(i.data(), i.node(), i.offset(), i.size());
    }
  delete[] data;
}

void
MPI2MemoryGrp::release_readwrite(void *data, distsize_t offset, int size)
{
  MemoryIter i(data, offsets_, n());
  for (i.begin(offset, size); i.ready(); i.next()) {
      replace_data(i.data(), i.node(), i.offset(), i.size());
    }
  delete[] data;
}

void
MPI2MemoryGrp::sync()
{
  MPI_Win_fence(0, rma_win_);
}

void *
MPI2MemoryGrp::localdata()
{
  return data_;
}

#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
