
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_memproc_h
#define _util_group_memproc_h

#include <stdio.h>
#include <sys/types.h>

#include <util/group/memmsg.h>

class ProcMemoryGrp: public MemoryGrp {
#define CLASSNAME ProcMemoryGrp
#define HAVE_KEYVAL_CTOR
#include <util/class/classd.h>
  private:
    char *data_;
  public:
    ProcMemoryGrp();
    ProcMemoryGrp(const RefKeyVal&);
    ~ProcMemoryGrp();

    void set_localsize(int);

    void *obtain_readwrite(int offset, int size);
    void *obtain_readonly(int offset, int size);
    void release_read(void *data, int offset, int size);
    void release_write(void *data, int offset, int size);

    void sync();
};

#endif
