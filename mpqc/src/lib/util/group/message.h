
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_message_h
#define _util_group_message_h

#include <util/state/state.h>
#include <util/state/classdRAVLMap.h>

class GrpReduce {
  public:
    virtual ~GrpReduce();
    virtual void reduce(void*target, void*data, int n, int size);
};

template <class T>
class GrpSumReduce: public GrpReduce {
  public:
    void reduce(void*target, void*data, int nelement, int element_size) {
        T *v1 = (T*)target, *v2 = (T*) data;
        for (int i=0; i<nelment; i++) {
            v1[i] += v2[i];
          }
      }
};

class MessageGrp {
  protected:
    // These are initialized by the initialize() member (see below).
    int me_;
    int n_;
    int nclass_;
    ClassDescPintRAVLMap classdesc_to_index_;
    ClassDescP *index_to_classdesc_;

    // The classdesc_to_index_ and index_to_classdesc_ arrays
    // cannot be initialized by the MessageGrp CTOR, because
    // the MessageGrp specialization has not yet been initialized
    // and communication is not available.  CTOR's of specializations
    // of MessageGrp must call the following member in their body
    // to complete the initialization process.
    void initialize(int me, int n);

    // Information about the last message received or probed.
    int last_source_;
    int last_size_; // the number of bytes
  public:
    MessageGrp();
    virtual ~MessageGrp();
    
    int me() { return me_; }
    int n() { return n_; }

    // send messages sequentially to the target
    virtual void send(int target, double* data, int ndata);
    virtual void send(int target, int* data, int ndata);
    virtual void send(int target, char* data, int nbyte);
    void send(int target, double data) { send(target,&data,1); }
    void send(int target, int data) { send(target,&data,1); }
    virtual void raw_send(int target, void* data, int nbyte) = 0;

    // send typed messages
    virtual void sendt(int target, int type, double* data, int ndata);
    virtual void sendt(int target, int type, int* data, int ndata);
    virtual void sendt(int target, int type, char* data, int nbyte);
    void sendt(int target, int type, double data) {sendt(target,type,&data,1);}
    void sendt(int target, int type, int data) {sendt(target,type&data,1);}
    virtual void raw_sendt(int target, int type, void* data, int nbyte) = 0;

    // receive message sent sequentually from the sender
    virtual void recv(int sender, double* data, int ndata);
    virtual void recv(int sender, int* data, int ndata);
    virtual void recv(int sender, char* data, int nbyte);
    void recv(int sender, double data) { recv(sender,&data,1); }
    void recv(int sender, int data) { recv(sender,&data,1); }
    virtual void raw_recv(int sender, void* data, int nbyte) = 0;

    // receive messages sent by type
    virtual void recvt(int type, double* data, int ndata);
    virtual void recvt(int type, int* data, int ndata);
    virtual void recvt(int type, char* data, int nbyte);
    void recvt(int type, double data) { recvt(type,&data,1); }
    void recvt(int type, int data) { recvt(type,&data,1); }
    virtual void raw_recvt(int type, void* data, int nbyte) = 0;

    // broadcast operations
    virtual void bcast(double* data, int ndata, int from = 0);
    virtual void bcast(int* data, int ndata, int from = 0);
    virtual void bcast(char* data, int nbyte, int from = 0);
    virtual void raw_bcast(void* data, int nbyte, int from);

    // global reduction operations
    virtual void reduce(double* data, int n, int target = -1);
    virtual void reduce(int* data, int n, int target = -1);
    virtual void reduce(char* data, int n, int target = -1);
    virtual void reduce(unsigned char* data, int n, int target = -1);
    virtual void raw_reduce(void*, int n, GrpReduce&, int target = -1);

    // synchronization
    virtual void sync() = 0;

    // Information about the last message received.
    int last_source() { return last_source_; }
    int last_size() { return last_size_; }

    // Each message group maintains an association of ClassDesc with
    // a global index so SavableState information can be sent between
    // nodes without needing to send the classname and look up the
    // ClassDesc with each transfer.  These routines return information
    // about that mapping.
    int classdesc_to_index(const ClassDesc*);
    const ClassDesc* index_to_classdesc(int);
};

// ProcMessageGrp provides a message group that has only one node
class ProcMessageGrp: public MessageGrp {
  public:
    ProcMessageGrp();
    ~ProcMessageGrp();
    void raw_send(int target, void* data, int nbyte);
    void raw_sendt(int target, int type, void* data, int nbyte);
    void raw_recv(int sender, void* data, int nbyte);
    void raw_recvt(int type, void* data, int nbyte);
    void raw_bcast(void* data, int nbyte, int from);
    void sync();
};

// The intMessageGrp uses integer message types to send and receive messages.
// The PICL and Paragon message groups derive from this.
class intMessageGrp: public MessageGrp {
  protected:
    int msgtype_nbit; // the total number of bits available
    const int ctl_nbit = 2; // control information bits
    int seq_nbit; // sequence information bits
    int typ_nbit; // type information bits
    int src_nbit; // source information bits

    // Masks for the fields in the type.
    const int ctl_mask = 0x3;
    int seq_mask;
    int typ_mask;
    int src_mask;

    // Shifts to construct a message type.
    int ctl_shift;
    int seq_shift;
    int typ_shift;
    int src_shift;

    int typ_msgtype(int source, int usrtype);
    int seq_msgtype(int source, int seq);

    // The next sequence number for each node is stored in these.
    int *source_seq;
    int *target_seq;
    
    void raw_send(int target, void* data, int nbyte);
    void raw_recv(int sender, void* data, int nbyte);
    void raw_sendt(int target, int type, void* data, int nbyte);
    void raw_recvt(int type, void* data, int nbyte);

    virtual void basic_send(int target, int type, void* data, int nbyte) = 0;
    virtual void basic_recv(int type, void* data, int nbyte) = 0;

    intMessageGrp();
    void initialize(int me, int n, int nbits);
  public:
    ~intMessageGrp();
};

class ShmMessageGrp: public intMessageGrp {
  protected:
    void basic_send(int target, int type, void* data, int nbyte);
    void basic_recv(int type, void* data, int nbyte);
    int basic_probe(int type);
    void initialize(int nprocs);
  public:
    ShmMessageGrp(); // read nprocs from environmental variable NUMPROC
    ShmMessageGrp(int nprocs);
    ~ShmMessageGrp();
    void sync();
};

#endif

