//
// messmpi.h
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

#ifndef _util_group_messmpi_h
#define _util_group_messmpi_h

#include <util/group/message.h>
#include <util/group/thread.h>

#define OMPI_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#include <mpi.h>

namespace sc {

/** The MPIMessageGrp class is an concrete implementation of MessageGrp
that uses the MPI 1 library.  */
class MPIMessageGrp: public MessageGrp {
  protected:
    void* buf;
    int bufsize;

    /// If true use the generic collective routines in the base class
    bool use_messagegrp_collectives_;

    /// Number of MPIMessageGrp's currently in use.
    static int nmpi_grps;
    /// lock to access nmpi_grps variable
    /// @note lifetime is managed manually to avoid destruction of lock before calling all MPIMessageGrp destructors
    static ThreadLock* grplock;
    /// Was MPI_Init called by one of MPIMessagrGrp? Will also call MPI_Finalize, if so
    static bool mpi_init_called;

    Ref<ThreadGrp> threadgrp;
    /// Currently each commgrp is a dup of MPI_COMM_WORLD
    MPI_Comm commgrp;
    
    /// Not thread-safe due to race condition on nmpi_grps variable.
    void init(MPI_Comm comm, int *argc=0, char ***argv=0);

    class MessageHandleData {
      public:
        MPI_Request req;
        size_t nbyte;
        MessageHandleData(size_t n): nbyte(n) {}
    };
  public:
    MPIMessageGrp();
    /** Use an MPI communicator to create a MessageGrp.  The comm
        argument could be a subset of MPI_COMM_WORLD, for example. */
    MPIMessageGrp(MPI_Comm comm);
    /** Use argc and argv to create a MPIMessageGrp.  This would
        have to be used for implementations of MPI that have MPI_Init
        fill in argc and argv. */
    MPIMessageGrp(int *argc, char ***argv);
    /** Construction MPIMessageGrp given a KeyVal input object. */
    MPIMessageGrp(const Ref<KeyVal>&);
    ~MPIMessageGrp();

    /// Clones (dups) an MPIMessageGrp from MPI_COMM_WORLD 
    Ref<MessageGrp> clone(void);
    Ref<MessageGrp> split(int grpkey=0, int rankkey=0);
    Ref<MessageGrp> subset(const std::set<int> &);
    
    void raw_send(int target, const void* data, int nbyte);
    void raw_recv(int sender, void* data, int nbyte,
                  MessageInfo *info=0);
    void raw_sendt(int target, int type, const void* data, int nbyte,
                   bool rcvrdy=false);
    void raw_recvt(int sender, int type, void* data, int nbyte,
                   MessageInfo *info=0);

    int probet(int sender, int type, MessageInfo *info=0);

    void sync();

    void sum(double*, int n, double*scratch = 0, int target = -1);
    void sum(int*, int n, int*scratch = 0, int target = -1);

    void reduce(double*, int n, GrpReduce<double>&,
                double*scratch = 0, int target = -1);
    void reduce(unsigned int*, int n, GrpReduce<unsigned int>&,
                unsigned int*scratch = 0, int target = -1);
    void reduce(int*, int n, GrpReduce<int>&,
                int*scratch = 0, int target = -1);
    void reduce(char*, int n, GrpReduce<char>&,
                char*scratch = 0, int target = -1);
    void reduce(unsigned char*, int n, GrpReduce<unsigned char>&,
                unsigned char*scratch = 0, int target = -1);
    void reduce(signed char*, int n, GrpReduce<signed char>&,
                signed char*scratch = 0, int target = -1);
    void reduce(short*, int n, GrpReduce<short>&,
                short*scratch = 0, int target = -1);
    void reduce(float*, int n, GrpReduce<float>&,
                float*scratch = 0, int target = -1);
    void reduce(long*, int n, GrpReduce<long>&,
                long*scratch = 0, int target = -1);

    void raw_nb_sendt(int sender, int type,
                      const void* data, int nbyte,
                      MessageHandle&,
                      bool rcvrdy=false);
    void raw_nb_recvt(int sender, int type,
                      void* data, int nbyte,
                      MessageHandle&);
    void wait(const MessageHandle&,
              MessageInfo *info=0);

    void raw_bcast(void* data, int nbyte, int from);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
