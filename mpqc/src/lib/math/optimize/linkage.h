
#ifndef _math_optimize_linkage_h
#define _math_optimize_linkage_h

#ifndef __PIC__

#include <math/optimize/qnewton.h>
#include <math/optimize/gdiis.h>
#include <math/optimize/efc.h>
#include <math/optimize/update.h>

#include <math/scmat/linkage.h>

const ClassDesc &optimize_force_link_a_ = QNewtonOpt::class_desc_;
const ClassDesc &optimize_force_link_b_ = GDIISOpt::class_desc_;
const ClassDesc &optimize_force_link_c_ = EFCOpt::class_desc_;
const ClassDesc &optimize_force_link_d_ = BFGSUpdate::class_desc_;
const ClassDesc &optimize_force_link_e_ = PowellUpdate::class_desc_;

#endif /* __PIC__ */

#endif
