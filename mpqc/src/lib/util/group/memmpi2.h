//
// memmpi2.h
//
// derived from code
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@sandia.gov>
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

#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_memmpi2_h
#define _util_group_memmpi2_h

#include <util/group/memamsg.h>

#include <mpi.h>

class MPI2MemoryGrp: public MsgMemoryGrp {
#define CLASSNAME MPI2MemoryGrp
#define HAVE_KEYVAL_CTOR
#include <util/class/classd.h>
  private:
    MPI_Win rma_win_;
    void *data_;
    
    void retrieve_data(void *, int node, int offset, int size);
    void replace_data(void *, int node, int offset, int size);
    void sum_data(double *data, int node, int doffset, int dsize);

  public:
    MPI2MemoryGrp(const RefMessageGrp& msg);
    MPI2MemoryGrp(const RefKeyVal&);
    ~MPI2MemoryGrp();

    void activate();
    void deactivate();

    void *localdata();
    void set_localsize(int);

    void *obtain_writeonly(distsize_t offset, int size);
    void *obtain_readwrite(distsize_t offset, int size);
    void *obtain_readonly(distsize_t offset, int size);
    void release_read(void *data, distsize_t offset, int size);
    void release_write(void *data, distsize_t offset, int size);
    void sum_reduction(double *data, distsize_t doffset, int dsize);
    void sum_reduction_on_node(double *data, int doffset, int dsize,
                               int node = -1);

    void sync();
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
