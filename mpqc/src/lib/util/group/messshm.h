
#ifndef _util_group_messshm_h
#define _util_group_messshm_h

#include <util/group/message.h>

//. \clsnm{ShmMessageGrp} is a concrete implementation of
//. \clsnmref{intMessageGrp}.
class ShmMessageGrp: public intMessageGrp {
#define CLASSNAME ShmMessageGrp
#define HAVE_KEYVAL_CTOR
#include <util/class/classd.h>
  protected:
    void basic_send(int target, int type, void* data, int nbyte);
    void basic_recv(int type, void* data, int nbyte);
    int basic_probe(int type);
    void initialize(int nprocs);
    void initialize();

    // Information about the last message received or probed.
    int last_type_;
    int last_source_;
    int last_size_; // the size in bytes

    void set_last_type(int a) { last_type_ = a; }
    void set_last_source(int a) { last_source_ = a; }
    void set_last_size(int a) { last_size_ = a; }
  public:
    ShmMessageGrp(); // read nprocs from environmental variable NUMPROC
    ShmMessageGrp(const RefKeyVal&);
    ShmMessageGrp(int nprocs);
    ~ShmMessageGrp();
    void sync();
 
    int last_source();
    int last_size();
    int last_type();
};

#endif
