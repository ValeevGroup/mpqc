
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

#include <chemistry/qc/dmtscf/scf_core.gbl>
#include <chemistry/qc/dmtscf/scf_inp.gbl>
#include <chemistry/qc/dmtscf/scf_dip.gbl>
#include <chemistry/qc/dmtscf/scf_esp.gbl>
#include <chemistry/qc/dmtscf/scf_lowd.gbl>
#include <chemistry/qc/dmtscf/scf_mull.gbl>
#include <chemistry/qc/dmtscf/scf_vect.gbl>
#include <chemistry/qc/dmtscf/scf_bnd.gbl>

#ifdef __cplusplus
}

/* The preferred input reader for C++ programs. */
int scf_get_keyval_input(KeyVal& keyval,
                         centers_t *_unique_centers,
                         scf_struct_t *_scf_info,
                         sym_struct_t *_sym_info,
                         scf_irreps_t *_irreps,
                         centers_t *_centers,
                         FILE *_outfile);
#endif /* __cplusplus */

#endif /* _scf_dmt_h */
