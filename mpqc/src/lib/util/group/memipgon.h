
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_memipgon_h
#define _util_group_memipgon_h

#include <stdio.h>
#include <util/group/memmid.h>

// This is a memory group that uses the paragon NX library,
// but doesn't use hrecv.  Instead irecv and isend are used.

class IParagonMemoryGrp: public MIDMemoryGrp {
#define CLASSNAME IParagonMemoryGrp
#include <util/class/classd.h>
  private:
    long lock();
    void unlock(long oldvalue);
    long send(void* data, int nbytes, int node, int type);
    long recv(void* data, int nbytes, int node, int type);
    long postrecv(void *data, int nbytes, int type);
    long wait(long, long = -1);
  public:
    IParagonMemoryGrp(const RefMessageGrp& msg, int localsize);
    ~IParagonMemoryGrp();
};

#endif
