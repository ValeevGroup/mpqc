
#ifndef _scf_dmt_h
#define _scf_dmt_

#ifdef __cplusplus
class KeyVal;
extern "C" {
#endif /* __cplusplus */

#include <chemistry/qc/dmtscf/scf.h>
#include <chemistry/qc/dmtscf/scfallc.h>
#include <chemistry/qc/dmtscf/scffree.h>
#include <chemistry/qc/dmtscf/scfzero.h>
#include <chemistry/qc/dmtscf/scfasgn.h>
#include <chemistry/qc/dmtscf/scfbrd.h>
#include <chemistry/qc/dmtscf/scfbwr.h>
#include <chemistry/qc/dmtscf/scfinit.h>
#include <chemistry/qc/dmtscf/scfiseq.h>
#include <chemistry/qc/dmtscf/scfprnt.h>
#include <chemistry/qc/dmtscf/scfbc0.h>
#include <chemistry/qc/dmtscf/scfrbc0.h>
#include <chemistry/qc/dmtscf/scfsbc0.h>
#include <chemistry/qc/dmtscf/scfrcv0.h>
#include <chemistry/qc/dmtscf/scfsnd0.h>

#include <chemistry/qc/dmtscf/scf_bnd.gbl>
#include <chemistry/qc/dmtscf/scf_core.gbl>
#include <chemistry/qc/dmtscf/scf_diis.gbl>
#include <chemistry/qc/dmtscf/scf_dip.gbl>
#include <chemistry/qc/dmtscf/scf_en.gbl>
#include <chemistry/qc/dmtscf/scf_esp.gbl>
#include <chemistry/qc/dmtscf/scf_ex.gbl>
#include <chemistry/qc/dmtscf/scf_fock.gbl>
#include <chemistry/qc/dmtscf/scf_gmat.gbl>
#include <chemistry/qc/dmtscf/scf_iter.gbl>
#include <chemistry/qc/dmtscf/scf_loopg.gbl>
#include <chemistry/qc/dmtscf/scf_loopj.gbl>
#include <chemistry/qc/dmtscf/scf_loopk.gbl>
#include <chemistry/qc/dmtscf/scf_lowd.gbl>
#include <chemistry/qc/dmtscf/scf_mkden.gbl>
#include <chemistry/qc/dmtscf/scf_mkgd.gbl>
#include <chemistry/qc/dmtscf/scf_mkgdlb.gbl>
#include <chemistry/qc/dmtscf/scf_mkj.gbl>
#include <chemistry/qc/dmtscf/scf_mkk.gbl>
#include <chemistry/qc/dmtscf/scf_mull.gbl>
#include <chemistry/qc/dmtscf/scf_oeis.gbl>
#include <chemistry/qc/dmtscf/scf_orth.gbl>
#include <chemistry/qc/dmtscf/scf_proj.gbl>
#include <chemistry/qc/dmtscf/scf_vect.gbl>

int scf_init_scf(centers_t*, scf_struct_t*, sym_struct_t*, char*);
int scf_init_scf_struct(centers_t*, scf_struct_t*, char*);

#ifdef __cplusplus
}

/* The preferred input reader for C++ programs. */
int scf_init_scf_struct_kv(KeyVal&, centers_t&, scf_struct_t&);
int scf_init_scf_kv(KeyVal&, centers_t&, scf_struct_t&, sym_struct_t&);

void scf_print_options(FILE*, scf_struct_t&);

#endif /* __cplusplus */

#endif /* _scf_dmt_h */



