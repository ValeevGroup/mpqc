
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_memmsg_h
#define _util_group_memmsg_h

#include <util/group/message.h>
#include <util/group/memory.h>

// A memory grp that initializes its data using a messagegrp.
class MsgMemoryGrp: public MemoryGrp {
#define CLASSNAME MsgMemoryGrp
#include <util/class/classda.h>

  protected:
    RefMessageGrp msg_;

  public:
    MsgMemoryGrp(const RefMessageGrp& msg, int localsize);
    ~MsgMemoryGrp();

    void sync();
};

#endif
