
#ifndef _math_scmat_linkage_h
#define _math_scmat_linkage_h

#ifndef __PIC__

#include <math/scmat/repl.h>
#include <math/scmat/dist.h>

#include <util/group/linkage.h>

const ClassDesc &math_scmat_force_link_a_ = ReplSCMatrixKit::class_desc_;
const ClassDesc &math_scmat_force_link_b_ = DistSCMatrixKit::class_desc_;

#endif /* __PIC__ */

#endif
