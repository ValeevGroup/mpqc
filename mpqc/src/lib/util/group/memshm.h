
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_memshm_h
#define _util_group_memshm_h

#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <iostream.h>

#include <util/group/globcnt.h>
#include <util/group/memmsg.h>

//. The \clsnm{ShmMessageGrp} class is an implementation of
//. \clsnmref{MessageGrp} that allows multiple process to be
//. started that communication with shared memory.
class ShmMemoryGrp: public MsgMemoryGrp {
#define CLASSNAME ShmMemoryGrp
#define HAVE_KEYVAL_CTOR
#include <util/class/classd.h>
  private:
    int shmid_;
    GlobalCounter lock_;
    GlobalCounter *update_;
    void *data_;
    void *memory_;
    Pool *pool_;
    RangeLock *rangelock_; // the locks_ member of the base class is ignored

    void clear_release_count();
    void wait_for_release();
    void note_release();
    void obtain_lock();
    void release_lock();

    void cleanup();
  public:
    ShmMemoryGrp(const RefMessageGrp& msg);
    ShmMemoryGrp(const RefKeyVal&);
    ~ShmMemoryGrp();

    void set_localsize(int);

    void *obtain_readwrite(int offset, int size);
    void *obtain_readonly(int offset, int size);
    void release_read(void *data, int offset, int size);
    void release_write(void *data, int offset, int size);

    virtual void sum_reduction(double *data, int doffset, int dsize);

    void print(ostream &o = cout);
};

#endif
