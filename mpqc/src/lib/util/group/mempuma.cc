//
// mempuma.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
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

#ifndef _util_group_mempuma_cc
#define _util_group_mempuma_cc

#ifdef __GNUC__
#pragma implementation
#endif

#include <util/group/mempuma.h>
#include <util/misc/formio.h>

extern "C" {
#include <nx.h>

int MPI2_RMA_init(void *base, int size, MPI_Datatype datatype, MPI_Comm comm,
                  MPI_Comm *newcomm );
int MPI2_Get(void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
             int target_rank, int target_displ, int target_count, 
             MPI_Datatype target_datatype, int target_increment,
             MPI_Comm comm );
 
int MPI2_Put(void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
             int target_rank, int target_displ, int target_count, 
             MPI_Datatype target_datatype, int target_increment,
             MPI_Comm comm );
}

#define DISABLE do { masktrap(1); cout.flush(); } while(0)
#define ENABLE do { cout.flush(); masktrap(0); } while(0)

#define PRINTF(args) do { DISABLE; \
                          cout << scprintf args; \
                          cout.flush(); \
                          ENABLE; \
                         } while(0)

#undef PRINTF
#define PRINTF(args)

///////////////////////////////////////////////////////////////////////
// The PumaMemoryGrp class

#define CLASSNAME PumaMemoryGrp
#define HAVE_KEYVAL_CTOR
#define PARENTS public ActiveMsgMemoryGrp
#include <util/class/classi.h>
void *
PumaMemoryGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  ActiveMsgMemoryGrp::_castdown(cd);
  return do_castdowns(casts,cd);
}

PumaMemoryGrp::PumaMemoryGrp(const RefMessageGrp& msg):
  ActiveMsgMemoryGrp(msg)
{
}

PumaMemoryGrp::PumaMemoryGrp(const RefKeyVal& keyval):
  ActiveMsgMemoryGrp(keyval)
{
}

PumaMemoryGrp::~PumaMemoryGrp()
{
  PRINTF(("%d: ~PumaMemoryGrp\n", me()));
}

void
PumaMemoryGrp::set_localsize(int localsize)
{
  if (debug_)
    cout << "PumaMemoryGrp::set_localsize(" << localsize << ")" << endl;

  ActiveMsgMemoryGrp::set_localsize(localsize);

  int r=MPI2_RMA_init(data_, localsize, MPI_BYTE, MPI_COMM_WORLD, &rma_comm_);
  if (r != MPI_SUCCESS) {
    cerr << scprintf("PumaMemoryGrp::set_localsize(%d) failed on %d",
                     localsize, me()) << endl;
    abort();
  }
}

void
PumaMemoryGrp::retrieve_data(void *data, int node, int offset, int size)
{
  PRINTF(("%d: retrieve_data: int offset = %d int size = %d\n",
          me(), offset/sizeof(int), size/sizeof(int)));

  int r = MPI2_Get(data, size, MPI_BYTE, node, offset, size, MPI_BYTE,
                   0, rma_comm_);
  if (r != MPI_SUCCESS) {
    cerr << "PumaMemoryGrp::retrieve_data failed" << endl;
    abort();
  }
}

void
PumaMemoryGrp::replace_data(void *data, int node, int offset, int size)
{
  PRINTF(("%d: replace_data: int offset = %d int size = %d\n",
          me(), offset/sizeof(int), size/sizeof(int)));

  int r = MPI2_Put(data, size, MPI_BYTE, node, offset, size, MPI_BYTE,
                   0, rma_comm_);
  if (r != MPI_SUCCESS) {
    cerr << "PumaMemoryGrp::replace_data failed" << endl;
    abort();
  }
}

void
PumaMemoryGrp::sum_data(double *data, int node, int offset, int size)
{
  int doffset = offset/sizeof(double);
  int dsize = size/sizeof(double);

  PRINTF(("%d: sum_data: doffset = %d dsize = %d node = %d\n",
          me(), doffset, dsize, node));

  double * tdata = new double[dsize];
  retrieve_data(tdata, node, offset, size);

  for (int i=0; i < dsize; i++) {
    tdata[i] += data[i];
  }

  replace_data(tdata,  node, offset, size);

  delete[] tdata;
}

#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
