
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_memshm_h
#define _util_group_memshm_h

#include <stdio.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>

#include <util/group/globcnt.h>
#include <util/group/memmsg.h>

#if defined(L486) || defined(PARAGON)
#ifndef SHMCTL_REQUIRES_SHMID
#  define SHMCTL_REQUIRES_SHMID
#endif
#endif

#if defined(L486) || defined(PARAGON)
#ifndef SHMDT_CHAR
#  define SHMDT_CHAR
#endif
#endif

class ShmMemoryGrp: public MsgMemoryGrp {
#define CLASSNAME ShmMemoryGrp
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
  public:
    ShmMemoryGrp(const RefMessageGrp& msg, int localsize);
    ~ShmMemoryGrp();

    void *obtain_readwrite(int offset, int size);
    void *obtain_readonly(int offset, int size);
    void release_read(void *data, int offset, int size);
    void release_write(void *data, int offset, int size);

    virtual void sum_reduction(double *data, int doffset, int dsize);

    void print(FILE *fp = stdout);
};

#endif
