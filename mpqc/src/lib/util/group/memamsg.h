
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_memamsg_h
#define _util_group_memamsg_h

#include <stdio.h>
#include <util/group/memmsg.h>

class MemoryLockRequest {
  public:
    enum { NData = 5 };
    enum Request { Deactivate, RelWrite, RelRead, RelReduce,
                   WriteOnly, ReadOnly, ReadWrite, Reduce };
  private:
    int data_[NData];
  public:
    MemoryLockRequest() {}
    MemoryLockRequest(Request r, int node = 0, int start = 0, int end = 0);
    void assign(Request r, int node, int start, int end);
    void *data() const { return (void *) data_; }
    int nbytes() const { return sizeof(int)*NData; }

    int request() const { return (Request) data_[0]; }
    int node() const { return data_[1]; }
    int start() const { return data_[2]; }
    int end() const { return data_[3]; }
    int serial_number() const { return data_[4]; }
};

class MemoryDataRequest {
  public:
    enum { NData = 5 };
    enum Request { Deactivate, Sync, Retrieve, Replace, DoubleSum };
  private:
    int data_[NData];
  public:
    MemoryDataRequest() {}
    MemoryDataRequest(Request r, int node = 0, int offset = 0, int size = 0);
    void assign(Request r, int node, int offset, int size);
    void *data() const { return (void *) data_; }
    int nbytes() const { return sizeof(int)*NData; }

    const char *request_string();

    MemoryDataRequest::Request request() const { return (Request) data_[0]; }
    int node() const { return data_[1]; }
    int offset() const { return data_[2]; }
    int size() const { return data_[3]; }
    int serial_number() const { return data_[4]; }

    void operator =(const MemoryDataRequest &r);

    void print(const char* msg = 0);
};

class MemoryDataRequestQueue {
  public:
    enum { MaxDepth = 20 };
  private:
    MemoryDataRequest q_[MaxDepth];
    int n_;
  public:
    MemoryDataRequestQueue(): n_(0) {}
    int n() const { return n_; }
    void push(MemoryDataRequest&);
    void pop(MemoryDataRequest&);

    MemoryDataRequest& operator[](int i) { return q_[i]; }
    void clear() { n_ = 0; }
};

class ActiveMsgMemoryGrp : public MsgMemoryGrp {
#define CLASSNAME ActiveMsgMemoryGrp
#include <util/class/classda.h>
  protected:
    char *data_;
    int use_locks_for_reduction_;

    // the defaults for these produce an error
    virtual void send_lock_request(MemoryLockRequest::Request,
                                   int offset, int size);
    virtual void wait_for_lock();

    virtual void retrieve_data(void *, int node, int offset, int size) = 0;
    virtual void replace_data(void *, int node, int offset, int size) = 0;
    virtual void sum_data(double *data, int node, int doffset, int dsize) = 0;
  public:
    ActiveMsgMemoryGrp(const RefMessageGrp& msg, int localsize);
    ~ActiveMsgMemoryGrp();

    void *obtain_writeonly(int offset, int size);
    void *obtain_readwrite(int offset, int size);
    void *obtain_readonly(int offset, int size);
    void release_read(void *data, int offset, int size);
    void release_write(void *data, int offset, int size);

    void sum_reduction(double *data, int doffset, int dsize);

    void print(FILE *fp = stdout);
};

#endif
