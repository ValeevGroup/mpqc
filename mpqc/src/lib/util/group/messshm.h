//
// messshm.h
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

#ifndef _util_group_messshm_h
#define _util_group_messshm_h

#include <unistd.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <sys/shm.h>

#include <util/group/message.h>

namespace sc {

#define SHMCOMMBUFSIZE 1500000

/* Set the maximum number of processors (including the host). */
#define MAXPROCS 17


struct commbuf_struct {
    int nmsg;
    int n_wait_for_change;
    int n_sync;
    char buf[SHMCOMMBUFSIZE];
};
typedef struct commbuf_struct commbuf_t;

struct msgbuf_struct {
    int type;
    int from;
    int size;
};
typedef struct msgbuf_struct msgbuf_t;

/** The ShmMessageGrp class is an implementation of MessageGrp that
allows multiple process to be started that communicate with shared memory.
This only provides improved performance if you have multiple CPU's in a
symmetric multiprocessor configuration.  Nonetheless, it is quite useful on
a single CPU for tracking down bugs.
*/
class ShmMessageGrp: public intMessageGrp {
  protected:
    void basic_send(int target, int type, const void* data, int nbyte);
    void basic_recv(int type, void* data, int nbyte);
    int basic_probe(int type);
    void initialize(int nprocs);
    void initialize();

    // previously static variables
    commbuf_t *commbuf[MAXPROCS];
    int shmid;
    int semid;
    int change_semid;
    void* sharedmem;
    struct sembuf semdec;
    struct sembuf seminc;

    // previously static functions for semephore operations
    msgbuf_t *NEXT_MESSAGE(msgbuf_t *m);
    void get_change(int node);
    void put_change(int node);
    void wait_for_write(int node);
    void release_write(int node);
#ifdef DEBUG
    void print_buffer(int node, int me);
#endif
  public:
    /// Reads the number of processors from environmental variable NUMPROC.
    ShmMessageGrp();
    /** The ShmMessageGrp KeyVal constructor takes a single keyword that
       specifies the number of processors.  Here is an example of a
       ParsedKeyVal input that creates a ShmMessageGrp that runs on four
       processors:

       <pre>
       message<ShmMessageGrp>: n = 4
       </pre>
    */
    ShmMessageGrp(const Ref<KeyVal>&);
    /// Initialize ShmMessageGrp to use nprocs processors.
    ShmMessageGrp(int nprocs);
    ~ShmMessageGrp();
    void sync();

    Ref<MessageGrp> clone(void);
};

}
     
#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
