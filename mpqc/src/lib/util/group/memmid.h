
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_memmid_h
#define _util_group_memmid_h

#include <util/group/memamsg.h>

// This is used for memory handler that use message identifiers to
// keep trace of messages.

class MIDMemoryGrp: public ActiveMsgMemoryGrp {
#define CLASSNAME MIDMemoryGrp
#include <util/class/classd.h>
  public:
    // This is public so memory handler functions can call it.
    void handler(long *mid = 0);
  protected:
    void handler(MemoryDataRequest&, long *mid = 0);
    
    void retrieve_data(void *, int node, int offset, int size);
    void replace_data(void *, int node, int offset, int size);
    void sum_data(double *data, int node, int doffset, int dsize);

    int active_;

    int data_request_type_;
    int data_type_to_handler_;
    int data_type_from_handler_;

    long data_request_mid_;

    MemoryDataRequest data_request_buffer_;

    int nsync_;

    int use_acknowledgments_;
    int use_active_messages_;

    void print_memreq(MemoryDataRequest &req,
                      const char * = 0, int target = -1);

    void do_wait(const char *msg, int mid,
                 MemoryDataRequestQueue &q, size_t expectedsize);
    void flush_queue(MemoryDataRequestQueue &q);

    virtual long lock() = 0;
    virtual void unlock(long oldvalue) = 0;
    virtual long send(void* data, int nbytes, int node, int type) = 0;
    virtual long recv(void* data, int nbytes, int node, int type) = 0;
    virtual long postrecv(void *data, int nbytes, int type) = 0;
    virtual long wait(long, long = -1) = 0;

    virtual void got_data_request_mid();
  public:
    MemoryDataRequest &data_request_buffer() { return data_request_buffer_; }

    MIDMemoryGrp(const RefMessageGrp& msg);
    MIDMemoryGrp(const RefKeyVal& msg);
    ~MIDMemoryGrp();

    void activate();
    void deactivate();

    void sync();
};

#endif
