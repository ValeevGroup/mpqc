
#ifndef _util_group_linkage_h
#define _util_group_linkage_h

#ifndef __PIC__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

# ifdef HAVE_SYSV_IPC
#   include <util/group/messshm.h>
    const ClassDesc &group_force_link_a_ = ShmMessageGrp::class_desc_;
# endif
const ClassDesc &fl9 = ProcMessageGrp::class_desc_;
# ifdef HAVE_NX_H
#  include <util/group/messpgon.h>
    const ClassDesc &group_force_link_b_ = ParagonMessageGrp::class_desc_;
# endif
#endif /* __PIC__ */

#endif
