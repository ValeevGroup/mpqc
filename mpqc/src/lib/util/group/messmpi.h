
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

#if HAVE_P4
    int nlocal;     // the number of processes on the master cluster
    int nremote;    // the number of remote clusters
    char *master;   // the name of the master cluster
    char * jobid;   // a unique job name selected by the user

    struct p4_cluster {
        char *hostname;  // name of the remote cluster
        int nslaves;     // the number of slaves on the remote cluster
    } *remote_clusters;

    struct p4_cluster * my_node_info(const char[], int&);
#endif
    
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
