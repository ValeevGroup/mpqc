
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_memmpl_h
#define _util_group_memmpl_h

#include <stdio.h>
#include <util/group/memamsg.h>

class MPLMemoryGrp: public MIDMemoryGrp {
#define CLASSNAME MPLMemoryGrp
#include <util/class/classd.h>
  private:
    long lock();
    void unlock(long oldvalue);
    long send(void* data, int nbytes, int node, int type);
    long recv(void* data, int nbytes, int node, int type);
    long postrecv(void *data, int nbytes, int type);
    long wait(long = -1);
  public:
    MPLMemoryGrp(const RefMessageGrp& msg, int localsize);
    ~MPLMemoryGrp();
};

#endif
