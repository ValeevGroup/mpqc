
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

GLOBAL_FUNCTION int
scf_oeis(_scf_info,_sym_info,_centers,S,T,V,H,_outfile)
scf_struct_t *_scf_info;
sym_struct_t *_sym_info;
centers_t *_centers;
dmt_matrix S;
dmt_matrix T;
dmt_matrix V;
dmt_matrix H;
FILE *_outfile;
{
  int i;
  int nlocal,li,lj;
  double *data;

/* initialize _centers struct */

  int_initialize_1e(0,0,_centers,_centers);
  int_initialize_offsets1(_centers,_centers);

/* calculate nuclear repulsion energy */

  _scf_info->nuc_rep = (double) int_nuclear_repulsion(_centers,_centers);
  if(mynode0()==0) {
    fprintf(_outfile,"\n  nuclear repulsion energy         = %f\n",
                                                     _scf_info->nuc_rep);
    fflush(_outfile);
    }


/* calculate one-electron integrals */

  nlocal = dmt_nlocal(S);

  for(i=0; i < nlocal ; i++) {
    dmt_get_block(S,i,&li,&lj,&data);
    int_shell_overlap(_centers,_centers,data,li,lj);
    }

  for(i=0; i < nlocal ; i++) {
    dmt_get_block(T,i,&li,&lj,&data);
    int_shell_kinetic(_centers,_centers,data,li,lj);
    }

  for(i=0; i < nlocal ; i++) {
    dmt_get_block(V,i,&li,&lj,&data);
    int_shell_nuclear(_centers,_centers,data,li,lj);
    }

  dmt_copy(V,H);
  dmt_sum(T,H);

 /* deallocate one-electron integral stuff */

  int_done_offsets1(_centers,_centers);
  int_done_1e();

  return(0);
  }
