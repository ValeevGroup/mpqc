
#ifndef _util_group_messshm_h
#define _util_group_messshm_h

#include <util/group/message.h>

class ParagonMessageGrp: public intMessageGrp {
#define CLASSNAME ParagonMessageGrp
#define HAVE_KEYVAL_CTOR
#include <util/class/classd.h>
  protected:
    void basic_send(int target, int type, void* data, int nbyte);
    void basic_recv(int type, void* data, int nbyte);
    int basic_probe(int type);
    void initialize();
  public:
    ParagonMessageGrp();
    ParagonMessageGrp(const RefKeyVal&);
    ~ParagonMessageGrp();
    void sync();
 
    int last_source();
    int last_size();
    int last_type();
};

#endif
