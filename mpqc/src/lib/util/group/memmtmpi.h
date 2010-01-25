//
// memmtmpi.h
// based on memmpi.h
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

#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_memmtmpi_h
#define _util_group_memmtmpi_h

#include <fstream>
#include <vector>

#define MPICH_SKIP_MPICXX
#include <mpi.h>

#include <util/group/message.h>
#include <util/group/memamsg.h>
#include <util/group/thread.h>
#include <util/misc/regtime.h>

namespace sc {

class MTMPIThread;

/** This MemoryGrp class requires a MT-safe MPI implementation.  The
default MessageGrp must be a MPIMessageGrp.  MPI must be safe with respect
to the default ThreadGrp.  Alternately, a MessageGrp and a ThreadGrp can be
passed to the constructor.  */
class MTMPIMemoryGrp: public ActiveMsgMemoryGrp {
  private:
    Ref<ThreadGrp> th_;

    Ref<ThreadLock> serial_lock_;
    int serial_;
    int serial(int node);

    MPI_Comm comp_comm_;
    MPI_Comm comm_comm_;
    int req_tag_;

    int active_;

    unsigned int *nreq_sent_;
    unsigned int *nreq_sent_buf_;

    MTMPIThread **thread_;
    Ref<ThreadLock> print_lock_; // needed for debugging only
    std::ofstream hout; // handler out
    std::ofstream mout; // main thread out

    void init_mtmpimg(MPI_Comm comm, int nthreads);

    Ref<RegionTimer> timer_;

    // Buffer data and manipulation members.
    int nbuffer_;

    int current_datareq_index_;
    std::vector<MemoryDataRequest> datareqs_;
    std::vector<MPI_Request> datareqs_mpireq_;

    int current_data_index_;
    std::vector<double*> databufs_;
    std::vector<MPI_Request> databufs_mpireq_;
    
    Ref<ThreadLock> buffer_lock_;

    int next_buffer(int &counter);
    void init_buffer();
    int next_buffer(int &counter,
                    std::vector<MPI_Request> &reqs);
    int get_buffer();
    int get_request();
    void done_buffers();

    // parent class pure virtuals
    void retrieve_data(void *, int node, long offset, long size, int lock);
    void replace_data(void *, int node, long offset, long size, int unlock);
    void sum_data(double *data, int node, long doffset, long dsize);

    friend class MTMPIThread;
  public:
    /** Construct a MTMPIMemoryGrp given a MessageGrp, ThreadGrp, and
        an MPI communicator.  The communicator can be a subset of
        MPI_COMM_WORLD, in which case, the MessageGrp must refer to the
        same subset. */
    MTMPIMemoryGrp(const Ref<MessageGrp>& msg, const Ref<ThreadGrp> &th,
                   MPI_Comm comm = MPI_COMM_WORLD);
    /** Construct a MTMPIMemoryGrp given a KeyVal input object. A
        fully thread safe MPI is needed (MPI_THREAD_MULTIPLE). The recognized
        keywords are:

        <table border="1">
        <tr><td>%Keyword<td>Type<td>Default<td>Description
        <tr><td><tt>num_threads</tt><td>integer<td>1<td>The number of threads to use for communication.
        <tr><td><tt>num_buffer</tt><td>integer<td>0<td>The number of buffers to prepost for communication.
        <tr><td><tt>use_timer</tt><td>boolean<td>false<td>Collect timing information.
        </table>
    */
    MTMPIMemoryGrp(const Ref<KeyVal> &);
    ~MTMPIMemoryGrp();

    void activate();
    void deactivate();

    void sync();

    Ref<MemoryGrp> clone();
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
