
#ifndef _chemistry_qc_scf_linkage_h
#define _chemistry_qc_scf_linkage_h

#ifndef __PIC__

#include <chemistry/qc/scf/clscf.h>
#include <chemistry/qc/scf/hsosscf.h>
#include <chemistry/qc/scf/ossscf.h>
#include <chemistry/qc/scf/tcscf.h>

#include <math/scmat/linkage.h>
#include <chemistry/molecule/linkage.h>

const ClassDesc &scf_force_link_a_ = CLSCF::class_desc_;
const ClassDesc &scf_force_link_b_ = HSOSSCF::class_desc_;
const ClassDesc &scf_force_link_c_ = OSSSCF::class_desc_;
const ClassDesc &scf_force_link_d_ = TCSCF::class_desc_;

#endif /* __PIC__ */

#endif
