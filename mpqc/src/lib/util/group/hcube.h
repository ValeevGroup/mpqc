
#ifndef _util_group_hcube_h
#define _util_group_hcube_h

#include <util/group/topology.h>

class HypercubeGMI: public GlobalMsgIter {
#define CLASSNAME HypercubeGMI
#include <util/class/classd.h>
  private:
    int nhyper_;
  protected:
    int fwdsendto();
    int fwdsend();
    int fwdrecvfrom();
    int fwdrecv();
  public:
    HypercubeGMI(int nproc, int me, int root);
    ~HypercubeGMI();
};

// This utilitizes a hypercube topology, but will work for any number of
// nodes.
class HypercubeTopology: public MachineTopology {
#define CLASSNAME HypercubeTopology
#define HAVE_KEYVAL_CTOR
#include <util/class/classd.h>
  public:
    HypercubeTopology();
    HypercubeTopology(const RefKeyVal&);
    ~HypercubeTopology();
    RefGlobalMsgIter global_msg_iter(const RefMessageGrp&, int target);
};

#endif
