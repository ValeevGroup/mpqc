
#ifndef _chemistry_qc_mbpt_mp2grad_h
#define _chemistry_qc_mbpt_mp2grad_h

#include <chemistry/qc/dmtscf/scf_dmt.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>
#include <math/array/math_lib.h>
#include <math/dmt/matrix.h>
#include <util/group/message.h>

void
mbpt_mp2_gradient(scf_struct_t &scf_info,
                  sym_struct_t &sym_info,
                  centers_t &centers,
                  int nfzc, int nfzv,
                  dmt_matrix Scf_Vec, dmt_matrix Fock, dmt_matrix FockO,
                  int mem_alloc,
                  FILE *outfile,
                  const RefMessageGrp &grp, const RefMemoryGrp &mem,
                  double &energy, double_matrix_t &gradient);

#endif
