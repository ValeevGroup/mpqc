
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_mempgon_h
#define _util_group_mempgon_h

#include <stdio.h>
#include <util/group/memamsg.h>

class ParagonMemoryGrp: public ActiveMsgMemoryGrp {
#define CLASSNAME ParagonMemoryGrp
#define HAVE_KEYVAL_CTOR
#include <util/class/classd.h>
  private:
    friend void paragon_memory_handler(long,long,long,long);
    
    void retrieve_data(void *, int node, int offset, int size);
    void replace_data(void *, int node, int offset, int size);
    void sum_data(double *data, int node, int doffset, int dsize);

    int data_request_type_;
    int data_type_to_handler_;
    int data_type_from_handler_;

    MemoryDataRequest data_request_buffer_;

    int active_;
  public:
    ParagonMemoryGrp(const RefMessageGrp& msg);
    ParagonMemoryGrp(const RefKeyVal&);
    ~ParagonMemoryGrp();

    void activate();
    void deactivate();
};

#endif
