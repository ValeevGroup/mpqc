
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_memmpl_h
#define _util_group_memmpl_h

#include <stdio.h>
#include <util/group/memamsg.h>

class MPLMemoryGrp: public ActiveMsgMemoryGrp {
#define CLASSNAME MPLMemoryGrp
#include <util/class/classd.h>
  private:
    friend static void mpl_memory_handler(int *mid);
    
    void retrieve_data(void *, int node, int offset, int size);
    void replace_data(void *, int node, int offset, int size);
    void sum_data(double *data, int node, int doffset, int dsize);

    int data_request_type_;
    int data_type_to_handler_;
    int data_type_from_handler_;

    MemoryDataRequest data_request_buffer_;

    int reactivate_;
    int active_;
    int nsync_;

    void print_memreq(MemoryDataRequest &req,
                      const char * = 0, int target = -1);
  public:
    MPLMemoryGrp(const RefMessageGrp& msg, int localsize);
    ~MPLMemoryGrp();

    void activate();
    void deactivate();

    void sync();
};

#endif
