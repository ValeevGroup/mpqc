//
// memmpi.h
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

#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_memmpi_h
#define _util_group_memmpi_h

#include <util/group/memmid.h>
#include <mpi.h>

class MPIMemoryGrp: public MIDMemoryGrp {
#define CLASSNAME MPIMemoryGrp
#define HAVE_KEYVAL_CTOR
#include <util/class/classd.h>
  private:
    enum { max_mid = 3 };
    int mid_ready_[max_mid];
    MPI_Request handles_[max_mid];

    long get_mid();
    void free_mid(long mid);
    void init_mid();

    long lockcomm();
    void unlockcomm(long oldvalue);
    long send(void* data, int nbytes, int node, int type);
    long recv(void* data, int nbytes, int node, int type);
    long postrecv(void *data, int nbytes, int type);
    long wait(long, long = -1);
    int probe(long);
  public:
    MPIMemoryGrp(const RefMessageGrp& msg);
    MPIMemoryGrp(const RefKeyVal &);
    ~MPIMemoryGrp();
    void deactivate();
};

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
