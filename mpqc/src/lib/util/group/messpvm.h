
#ifndef _util_group_messpvm_h
#define _util_group_messpvm_h

#include <util/group/message.h>

class PVMMessageGrp: public MessageGrp {
#define CLASSNAME PVMMessageGrp
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classd.h>
  protected:
    int *tids;

    int rtid;
    int rtag;
    int rlen;
  public:
    PVMMessageGrp();
    PVMMessageGrp(const RefKeyVal&);
    ~PVMMessageGrp();

    void raw_send(int target, void* data, int nbyte);
    void raw_recv(int sender, void* data, int nbyte);
    void raw_sendt(int target, int type, void* data, int nbyte);
    void raw_recvt(int type, void* data, int nbyte);

    int probet(int type);

    int last_source();
    int last_size();
    int last_type();
};

#endif
