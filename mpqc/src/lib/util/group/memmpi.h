
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_memmpi_h
#define _util_group_memmpi_h

#include <stdio.h>
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

    long lock();
    void unlock(long oldvalue);
    long send(void* data, int nbytes, int node, int type);
    long recv(void* data, int nbytes, int node, int type);
    long postrecv(void *data, int nbytes, int type);
    long wait(long, long = -1);
  public:
    MPIMemoryGrp(const RefMessageGrp& msg);
    MPIMemoryGrp(const RefKeyVal &);
    ~MPIMemoryGrp();
    void deactivate();
};

#endif
