
#include <stdio.h>
#include <tmpl.h>
#include <comm/picl/picl.h>
#include <comm/picl/ext/piclext.h>
#include <math/dmt/libdmt.h>
#include <math/array/math_lib.h>
#include <util/bio/libbio.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>
#include "scf.h"

#include "scf_oeis.gbl"
#include "scf_oeis.lcl"

/*************************************************************************
 *
 * given a centers struct and four scattered matrices, calculate the
 * one electron integral
 *
 * input:
 *   scf_info = pointer to initialized scf struct
 *   centers  = pointer to initialized centers struct
 *   S, T, V, H = scattered dmt matrices
 *   outfile  = FILE pointer to output (may be null)
 *
 * on return:
 *   scf_info contains nuclear repulsion energy
 *   S contains overlap integrals
 *   T contains kinetic energy integrals
 *   V contains potential energy (nucl-elect attrac) integrals
 *   H contains core hamiltonian (T+V)
 *
 * return 0 on success, -1 on failure
 */

GLOBAL_FUNCTION int
scf_oeis(scf_info, centers, S, T, V, H, outfile)
scf_struct_t *scf_info;
centers_t *centers;
dmt_matrix S;
dmt_matrix T;
dmt_matrix V;
dmt_matrix H;
FILE *outfile;
{
  int i;
  int nlocal,li,lj;
  double *data;

/* initialize _centers struct */

  int_initialize_1e(0,0,centers,centers);
  int_initialize_offsets1(centers,centers);

/* calculate nuclear repulsion energy */

  scf_info->nuc_rep = (double) int_nuclear_repulsion(centers,centers);
  if (mynode0()==0 && outfile) {
    fprintf(outfile,"\n  nuclear repulsion energy         = %f\n",
            scf_info->nuc_rep);
    fflush(outfile);
  }


 /* calculate one-electron integrals */
  nlocal = dmt_nlocal(S);

  for (i=0; i < nlocal ; i++) {
    dmt_get_block(S,i,&li,&lj,&data);
    int_shell_overlap(centers,centers,data,li,lj);
  }

  for (i=0; i < nlocal ; i++) {
    dmt_get_block(T,i,&li,&lj,&data);
    int_shell_kinetic(centers,centers,data,li,lj);
  }

  for (i=0; i < nlocal ; i++) {
    dmt_get_block(V,i,&li,&lj,&data);
    int_shell_nuclear(centers,centers,data,li,lj);
  }

  dmt_copy(V,H);
  dmt_sum(T,H);

 /* deallocate one-electron integral stuff */

  int_done_offsets1(centers,centers);
  int_done_1e();

  return(0);
}
