
#ifndef _util_group_topology_h
#define _util_group_topology_h

#include <util/class/class.h>
#include <util/keyval/keyval.h>

DescribedClass_REF_fwddec(MessageGrp);
class GlobalMsgIter: public DescribedClass {
#define CLASSNAME GlobalMsgIter
#include <util/class/classda.h>
  protected:
    int me_;
    int nproc_;
    int root_;
    int i_;
    int n_; // the number of steps--intialized by derived class CTORs
    int fwd_;

    // for sending messages in the forward direction (like a bcast)
    virtual int fwdsendto() = 0;
    virtual int fwdsend() = 0;
    virtual int fwdrecvfrom() = 0;
    virtual int fwdrecv() = 0;
  public:
    GlobalMsgIter(int nproc, int me, int root = 0);
    void backwards() { fwd_ = 0; i_ = n_ - 1; }
    void forwards() { fwd_ = 1; i_ = 0; }
    void next() { if (fwd_) i_++; else i_--; }
    int done() { return i_<0 || i_>=n_; }
    int n() { return n_; }
    int sendto() { return fwd_?fwdsendto():fwdrecvfrom(); }
    int send() { return fwd_?fwdsend():fwdrecv(); }
    int recvfrom() { return fwd_?fwdrecvfrom():fwdsendto(); }
    int recv() { return fwd_?fwdrecv():fwdsend(); }
};
DescribedClass_REF_dec(GlobalMsgIter);

class MachineTopology: public DescribedClass {
#define CLASSNAME MachineTopology
#include <util/class/classda.h>
  public:
    MachineTopology();
    MachineTopology(const RefKeyVal&);
    ~MachineTopology();

    virtual RefGlobalMsgIter global_msg_iter(const RefMessageGrp&,
                                             int target) = 0;
};
DescribedClass_REF_dec(MachineTopology);

#endif
