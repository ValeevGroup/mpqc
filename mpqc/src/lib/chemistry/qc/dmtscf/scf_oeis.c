
/* $Log$
 * Revision 1.1  1993/12/29 12:53:16  etseidl
 * Initial revision
 *
 * Revision 1.5  1992/06/17  21:54:27  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.4  1992/05/13  18:23:22  jannsen
 * Added fflush calls.
 *
 * Revision 1.3  1992/05/12  10:35:48  seidl
 * calculate nuclear rep energy here now
 *
 * Revision 1.2  1992/03/21  00:40:54  seidl
 * change sym_libv2.h to chemistry/qc/dmtsym/sym_dmt.h
 *
 * Revision 1.1.1.1  1992/03/17  16:26:50  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:26:49  seidl
 * Initial revision
 *
 * Revision 1.2  1992/02/18  17:50:43  seidl
 * no longer allocate oei matrices here
 *
 * Revision 1.1  1992/02/04  23:48:08  seidl
 * Initial revision
 *
 * Revision 1.5  1992/01/16  19:53:07  seidl
 * use new integral routines
 *
 * Revision 1.4  1992/01/09  11:48:38  seidl
 * add parallel code, master file pointers initialized in scf_open_files
 *
 * Revision 1.3  1991/12/20  16:31:11  seidl
 * fix bug in sq_to_tri
 *
 * Revision 1.2  1991/12/17  21:44:48  seidl
 * change to new util/sgen/sgen format
 *
 * Revision 1.1  1991/12/02  19:58:51  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#include <stdio.h>
#include <tmpl.h>
#include <math/dmt/libdmt.h>
#include <math/array/math_lib.h>
#include <util/bio/libbio.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>
#include "scf.h"

#include "scf_oeis.gbl"
#include "scf_oeis.lcl"

GLOBAL_FUNCTION int
scf_oeis(_scf_info,_sym_info,_irreps,_centers,S,T,V,H,_outfile)
scf_struct_t *_scf_info;
sym_struct_t *_sym_info;
scf_irreps_t *_irreps;
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
