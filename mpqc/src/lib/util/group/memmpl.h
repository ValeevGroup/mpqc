
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_memmpl_h
#define _util_group_memmpl_h

#include <util/group/memmid.h>

class MPLMemoryGrp: public MIDMemoryGrp {
#define CLASSNAME MPLMemoryGrp
#define HAVE_KEYVAL_CTOR
#include <util/class/classd.h>
  private:
    long lockcomm();
    void unlockcomm(long oldvalue);
    long send(void* data, int nbytes, int node, int type);
    long recv(void* data, int nbytes, int node, int type);
    long postrecv(void *data, int nbytes, int type);
    long wait(long, long = -1);
  public:
    MPLMemoryGrp(const RefMessageGrp& msg);
    MPLMemoryGrp(const RefKeyVal&);
    ~MPLMemoryGrp();
    void deactivate();
};

#endif
