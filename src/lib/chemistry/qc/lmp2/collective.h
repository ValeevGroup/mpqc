
/*
 * Copyright 2009 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 *
 * This file is a part of the MPQC LMP2 library.
 *
 * The MPQC LMP2 library is free software: you can redistribute it
 * and/or modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _chemistry_qc_lmp2_collective_h
#define _chemistry_qc_lmp2_collective_h

#include <mpqc_config.h>

#ifdef HAVE_MPI

#define OMPI_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#include <mpi.h>

namespace sc {

/// Implements an all-to-all communication operation.
/// This may be in terms of MPI_Alltoallv or a custom
/// implementation that is often faster.
void custom_alltoallv(void *sendbuf,
                      int *sendcnts, 
                      int *sdispls, 
                      MPI_Datatype sendtype,
                      void *recvbuf, 
                      int *recvcnts, 
                      int *rdispls,
                      MPI_Datatype recvtype,
                      MPI_Comm comm);

/// Implements a pipelined reduce-broadcast operation.
/// This is often faster than that provided by MPI_Allreduce.
void allsum(double *data, long n);

}

#endif


#endif

