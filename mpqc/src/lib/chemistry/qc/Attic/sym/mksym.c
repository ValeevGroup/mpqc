
/* $Log$
 * Revision 1.1  1993/12/29 12:53:13  etseidl
 * Initial revision
 *
 * Revision 1.1.1.1  1992/03/17  17:10:57  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:10:55  seidl
 * Initial revision
 *
 * Revision 1.2  1992/01/21  11:32:32  seidl
 * use libintv2
 *
 * Revision 1.1  1991/12/20  15:45:23  seidl
 * Initial revision
 *
 * Revision 1.1  1991/12/03  16:42:20  etseidl
 * Initial revision
 *
 * Revision 1.1  1991/11/22  18:28:57  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#include <stdio.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <util/ipv2/ip_libv2.h>
#include <chemistry/qc/intv2/int_libv2.h>

#include "symm.h"
#include "create_r.gbl"

#include "mksym.gbl"
#include "mksym.lcl"

GLOBAL_FUNCTION VOID
sym_sym_matrix(_centers,_sym_info,skel,sym,tr,scr,_outfile)
centers_t *_centers;
sym_struct_t *_sym_info;
double_matrix_t *skel;
double_matrix_t *sym;
double_matrix_t *tr;
double_matrix_t *scr;
FILE *_outfile;
{
  int i,j;
  int ng;

  for(i=0; i < sym->n1; i++) bzero(sym->d[i],sizeof(double)*sym->n2);

/* Remember, tr->d holds the transpose of the R matrix defined in
 * Dupuis & King's paper */

  for(ng=0; ng < _sym_info->g ; ng++) {
    sym_create_r(_centers,_sym_info,tr,ng,_outfile);

    math_dmxdm_dm(skel,0,tr,1,scr,1,skel->n1,skel->n1,skel->n1,0);
    math_dmxdm_dm(tr,0,scr,1,sym,0,tr->n1,tr->n1,tr->n1,1);
    }

  for(i=0; i < _centers->nfunc ; i++)
    for(j=0; j < _centers->nfunc ; j++)
      sym->d[i][j] /= _sym_info->g;
  }

GLOBAL_FUNCTION VOID
sym_sym_vector(_centers,_sym_info,skel,sym,tr,scr,_outfile)
centers_t *_centers;
sym_struct_t *_sym_info;
double_vector_t *skel;
double_vector_t *sym;
double_matrix_t *tr;
double_matrix_t *scr;
FILE *_outfile;
{
  int i,j;
  int ng;

  bzero(sym->d,sizeof(double)*sym->n);

/* Remember, tr->d holds the transpose of the R matrix defined in
 * Dupuis & King's paper */

  for(ng=0; ng < _sym_info->g ; ng++) {
    sym_create_r(_centers,_sym_info,tr,ng,_outfile);

    math_dvxdm_dm(skel,0,tr,1,scr,1,tr->n1,tr->n1,tr->n1,0);
    math_dmxdm_dv(tr,0,scr,1,sym,0,tr->n1,tr->n1,tr->n1,1);
    }

  for(i=0; i < sym->n; i++) sym->d[i] /= _sym_info->g;
  }
