
#ifndef _util_group_messmpi_h
#define _util_group_messmpi_h

#include <util/group/message.h>

class MPIMessageGrp: public MessageGrp {
#define CLASSNAME MPIMessageGrp
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classd.h>
  protected:
    void* buf;
    int bufsize;

    int rnode;
    int rtag;
    int rlen;

    void init();
  public:
    MPIMessageGrp();
    MPIMessageGrp(const RefKeyVal&);
    ~MPIMessageGrp();

    void raw_send(int target, void* data, int nbyte);
    void raw_recv(int sender, void* data, int nbyte);
    void raw_sendt(int target, int type, void* data, int nbyte);
    void raw_recvt(int type, void* data, int nbyte);

    int probet(int type);

    int last_source();
    int last_size();
    int last_type();

    void sync();

    void reduce(double*, int n, GrpReduce<double>&,
                double*scratch = 0, int target = -1);

    void raw_bcast(void* data, int nbyte, int from);
};

#endif
