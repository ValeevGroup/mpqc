
#ifndef _chemistry_qc_mbpt_ffo_h
#define _chemistry_qc_mbpt_ffo_h

#include <math/dmt/matrix.h>
#include <chemistry/qc/dmtscf/scf_dmt.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>

int
mbpt_ffo(dmt_matrix S,
         scf_struct_t *_scf_info, sym_struct_t *_sym_info, centers_t *_centers,
         dmt_matrix SCF_VEC, dmt_matrix FOCK, dmt_matrix FOCKO);

#endif
