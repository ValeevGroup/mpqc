
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_message_h
#define _util_group_message_h

#include <math.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/state/classdRAVLMap.h>
#include <util/keyval/keyval.h>
#include <util/group/topology.h>

template <class T>
class GrpReduce {
  public:
    virtual ~GrpReduce() {};
    virtual void reduce(T*target, T*data, int n) = 0;
};

template <class T>
class GrpSumReduce: public GrpReduce<T> {
  public:
    void reduce(T*target, T*data, int nelement) {
        for (int i=0; i<nelement; i++) {
            target[i] += data[i];
          }
      }
};

template <class T>
class GrpMinReduce: public GrpReduce<T> {
  public:
    void reduce(T*target, T*data, int nelement) {
        for (int i=0; i<nelement; i++) {
            if (target[i] > data[i]) target[i] = data[i];
          }
      }
};

template <class T>
class GrpMaxReduce: public GrpReduce<T> {
  public:
    void reduce(T*target, T*data, int nelement) {
        for (int i=0; i<nelement; i++) {
            if (target[i] < data[i]) target[i] = data[i];
          }
      }
};

template <class T>
class GrpArithmeticAndReduce: public GrpReduce<T> {
  public:
    void reduce(T*target, T*data, int nelement) {
        for (int i=0; i<nelement; i++) {
            target[i] = target[i] & data[i];
          }
      }
};

template <class T>
class GrpArithmeticOrReduce: public GrpReduce<T> {
  public:
    void reduce(T*target, T*data, int nelement) {
        for (int i=0; i<nelement; i++) {
            target[i] = target[i] | data[i];
          }
      }
};

template <class T>
class GrpArithmeticXOrReduce: public GrpReduce<T> {
  public:
    void reduce(T*target, T*data, int nelement) {
        for (int i=0; i<nelement; i++) {
            target[i] = target[i] ^ data[i];
          }
      }
};

template <class T>
class GrpProductReduce: public GrpReduce<T> {
  public:
    void reduce(T*target, T*data, int nelement) {
        for (int i=0; i<nelement; i++) {
            target[i] *= data[i];
          }
      }
};

template <class T>
class GrpFunctionReduce: public GrpReduce<T> {
  private:
    void (*func_)(T*target,T*data,int nelement);
  public:
    GrpFunctionReduce(void(*func)(T*,T*,int)):func_(func) {}
    void reduce(T*target, T*data, int nelement) {
        (*func_)(target,data,nelement);
      }
};

DescribedClass_REF_fwddec(MessageGrp);

class MessageGrp: public DescribedClass {
#define CLASSNAME MessageGrp
#include <util/class/classda.h>
  private:
    // These are initialized by the initialize() member (see below).
    int me_;
    int n_;
    int nclass_;
    ClassDescPintRAVLMap classdesc_to_index_;
    ClassDescP *index_to_classdesc_;
  protected:
    // The classdesc_to_index_ and index_to_classdesc_ arrays
    // cannot be initialized by the MessageGrp CTOR, because
    // the MessageGrp specialization has not yet been initialized
    // and communication is not available.  CTOR's of specializations
    // of MessageGrp must call the following member in their body
    // to complete the initialization process.
    void initialize(int me, int n);

    RefMachineTopology topology_;
  public:
    MessageGrp();
    MessageGrp(const RefKeyVal&);
    virtual ~MessageGrp();
    
    int me() { return me_; }
    int n() { return n_; }

    // The default message group contains the primary message group to
    // be used by an application.
    static void set_default_messagegrp(const RefMessageGrp&);
    static MessageGrp* get_default_messagegrp();

    // The initial message group is the group that starts up a process.
    // This returns null if this process is first and it is up to the
    // programmer to create a messagegrp.
    static MessageGrp* initial_messagegrp();

    // send messages sequentially to the target
    virtual void send(int target, double* data, int ndata);
    virtual void send(int target, int* data, int ndata);
    virtual void send(int target, char* data, int nbyte);
    virtual void send(int target, unsigned char* data, int nbyte);
    virtual void send(int target, short* data, int ndata);
    virtual void send(int target, long* data, int ndata);
    virtual void send(int target, float* data, int ndata);
    void send(int target, double data) { send(target,&data,1); }
    void send(int target, int data) { send(target,&data,1); }
    virtual void raw_send(int target, void* data, int nbyte) = 0;

    // send typed messages
    virtual void sendt(int target, int type, double* data, int ndata);
    virtual void sendt(int target, int type, int* data, int ndata);
    virtual void sendt(int target, int type, char* data, int nbyte);
    virtual void sendt(int target, int type, unsigned char* data, int nbyte);
    virtual void sendt(int target, int type, short* data, int ndata);
    virtual void sendt(int target, int type, long* data, int ndata);
    virtual void sendt(int target, int type, float* data, int ndata);
    void sendt(int target, int type, double data) {sendt(target,type,&data,1);}
    void sendt(int target, int type, int data) {sendt(target,type&data,1);}
    virtual void raw_sendt(int target, int type, void* data, int nbyte) = 0;

    // receive message sent sequentually from the sender
    virtual void recv(int sender, double* data, int ndata);
    virtual void recv(int sender, int* data, int ndata);
    virtual void recv(int sender, char* data, int nbyte);
    virtual void recv(int sender, unsigned char* data, int nbyte);
    virtual void recv(int sender, short* data, int ndata);
    virtual void recv(int sender, long* data, int ndata);
    virtual void recv(int sender, float* data, int ndata);
    void recv(int sender, double data) { recv(sender,&data,1); }
    void recv(int sender, int data) { recv(sender,&data,1); }
    virtual void raw_recv(int sender, void* data, int nbyte) = 0;

    // receive messages sent by type
    virtual void recvt(int type, double* data, int ndata);
    virtual void recvt(int type, int* data, int ndata);
    virtual void recvt(int type, char* data, int nbyte);
    virtual void recvt(int type, unsigned char* data, int nbyte);
    virtual void recvt(int type, short* data, int ndata);
    virtual void recvt(int type, long* data, int ndata);
    virtual void recvt(int type, float* data, int ndata);
    void recvt(int type, double data) { recvt(type,&data,1); }
    void recvt(int type, int data) { recvt(type,&data,1); }
    virtual void raw_recvt(int type, void* data, int nbyte) = 0;

    // ask if a given type message has been received
    virtual int probet(int type) = 0;

    // broadcast operations
    virtual void bcast(double* data, int ndata, int from = 0);
    virtual void bcast(int* data, int ndata, int from = 0);
    virtual void bcast(char* data, int nbyte, int from = 0);
    virtual void bcast(unsigned char* data, int nbyte, int from = 0);
    virtual void bcast(short* data, int ndata, int from = 0);
    virtual void bcast(long* data, int ndata, int from = 0);
    virtual void bcast(float* data, int ndata, int from = 0);
    virtual void raw_bcast(void* data, int nbyte, int from);

    // global reduction operations
    virtual void sum(double* data, int n, double* = 0, int target = -1);
    virtual void sum(int* data, int n, int* = 0, int target = -1);
    virtual void sum(char* data, int n, char* = 0, int target = -1);
    virtual void sum(unsigned char* data, int n,
                     unsigned char* = 0, int target = -1);
    virtual void max(double* data, int n, double* = 0, int target = -1);
    virtual void max(int* data, int n, int* = 0, int target = -1);
    virtual void max(char* data, int n, char* = 0, int target = -1);
    virtual void max(unsigned char* data, int n,
                     unsigned char* = 0, int target = -1);
    virtual void min(double* data, int n, double* = 0, int target = -1);
    virtual void min(int* data, int n, int* = 0, int target = -1);
    virtual void min(char* data, int n, char* = 0, int target = -1);
    virtual void min(unsigned char* data, int n,
                     unsigned char* = 0, int target = -1);
    virtual void reduce(double*, int n, GrpReduce<double>&,
                        double*scratch = 0, int target = -1);
    virtual void reduce(int*, int n, GrpReduce<int>&,
                        int*scratch = 0, int target = -1);
    virtual void reduce(char*, int n, GrpReduce<char>&,
                        char*scratch = 0, int target = -1);
    virtual void reduce(unsigned char*, int n, GrpReduce<unsigned char>&,
                        unsigned char*scratch = 0, int target = -1);
    virtual void reduce(short*, int n, GrpReduce<short>&,
                        short*scratch = 0, int target = -1);
    virtual void reduce(float*, int n, GrpReduce<float>&,
                        float*scratch = 0, int target = -1);
    virtual void reduce(long*, int n, GrpReduce<long>&,
                        long*scratch = 0, int target = -1);

    // synchronization
    virtual void sync();

    // Information about the last message received or probed.
    virtual int last_source() = 0;
    virtual int last_size() = 0;
    virtual int last_type() = 0;

    // Each message group maintains an association of ClassDesc with
    // a global index so SavableState information can be sent between
    // nodes without needing to send the classname and look up the
    // ClassDesc with each transfer.  These routines return information
    // about that mapping.
    int classdesc_to_index(const ClassDesc*);
    const ClassDesc* index_to_classdesc(int);
};
DescribedClass_REF_dec(MessageGrp);

// ProcMessageGrp provides a message group that has only one node
class ProcMessageGrp: public MessageGrp {
#define CLASSNAME ProcMessageGrp
#define HAVE_KEYVAL_CTOR
#include <util/class/classd.h>
  private:
    // Information about the last message received or probed.
    int last_type_;
    int last_source_;
    int last_size_; // the size in bytes

    void set_last_type(int a) { last_type_ = a; }
    void set_last_source(int a) { last_source_ = a; }
    void set_last_size(int a) { last_size_ = a; }
  public:
    ProcMessageGrp();
    ProcMessageGrp(const RefKeyVal&);
    ~ProcMessageGrp();
    void raw_send(int target, void* data, int nbyte);
    void raw_sendt(int target, int type, void* data, int nbyte);
    void raw_recv(int sender, void* data, int nbyte);
    void raw_recvt(int type, void* data, int nbyte);
    void raw_bcast(void* data, int nbyte, int from);
    int probet(int type);
    void sync();
 
    int last_source();
    int last_size();
    int last_type();
};

// The intMessageGrp uses integer message types to send and receive messages.
// The PICL and Paragon message groups derive from this.
class intMessageGrp: public MessageGrp {
#define CLASSNAME intMessageGrp
#include <util/class/classda.h>
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

    int msgtype_typ(int msgtype);
    int typ_msgtype(int usrtype);
    int seq_msgtype(int source, int seq);

    // The next sequence number for each node is stored in these.
    int *source_seq;
    int *target_seq;
    
    virtual void basic_send(int target, int type, void* data, int nbyte) = 0;
    virtual void basic_recv(int type, void* data, int nbyte) = 0;
    virtual int basic_probe(int type) = 0;

    intMessageGrp();
    intMessageGrp(const RefKeyVal&);
    void initialize(int me, int n, int nbits);
  public:
    ~intMessageGrp();

    void raw_send(int target, void* data, int nbyte);
    void raw_recv(int sender, void* data, int nbyte);
    void raw_sendt(int target, int type, void* data, int nbyte);
    void raw_recvt(int type, void* data, int nbyte);

    int probet(int);
};

#endif

