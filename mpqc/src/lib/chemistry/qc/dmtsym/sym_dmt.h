#ifndef _sym_dmt_h
#define _sym_dmt_h

#ifdef __cplusplus

extern "C" {
#endif

#include <math/dmt/libdmt.h>

#include <chemistry/qc/intv2/atoms.h>

#include <chemistry/qc/dmtsym/symm.h>
#include <chemistry/qc/dmtsym/symm_mac.h>
#include <chemistry/qc/dmtsym/symminit.h>
#include <chemistry/qc/dmtsym/symmallc.h>
#include <chemistry/qc/dmtsym/symmasgn.h>
#include <chemistry/qc/dmtsym/symmiseq.h>
#include <chemistry/qc/dmtsym/symmzero.h>
#include <chemistry/qc/dmtsym/symmfree.h>
#include <chemistry/qc/dmtsym/symmprnt.h>
#include <chemistry/qc/dmtsym/symmbc0.h>
#include <chemistry/qc/dmtsym/symmrbc0.h>
#include <chemistry/qc/dmtsym/symmsbc0.h>
#include <chemistry/qc/dmtsym/symmrcv0.h>
#include <chemistry/qc/dmtsym/symmsnd0.h>

#include <chemistry/qc/dmtsym/create_r.gbl>
#include <chemistry/qc/dmtsym/mksym.gbl>
#include <chemistry/qc/dmtsym/skeleton.gbl>
#include <chemistry/qc/dmtsym/syminit.gbl>

#ifdef __cplusplus
}

#include <util/keyval/keyval.h>
#include <math/symmetry/pointgrp.h>
#include <chemistry/qc/basis/gaussbas.h>

int sym_struct_from_gbs(const RefGaussianBasisSet&, sym_struct_t&);
int sym_struct_from_pg(const PointGroup&, centers_t&, sym_struct_t&);
int sym_init_centers(const RefKeyVal&, centers_t&, sym_struct_t&);

#endif

#endif
