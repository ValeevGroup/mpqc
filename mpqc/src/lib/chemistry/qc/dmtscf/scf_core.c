
/* $Log$
 * Revision 1.1  1993/12/29 12:53:14  etseidl
 * Initial revision
 *
 * Revision 1.6  1992/06/29  17:47:11  seidl
 * change scr to upper case
 *
 * Revision 1.5  1992/06/24  15:25:38  seidl
 * change dmt_matrices to upper case
 *
 * Revision 1.4  1992/06/17  21:54:05  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.3  1992/04/22  15:53:35  seidl
 * do not free sahalf if proj_vector is true
 *
 * Revision 1.2  1992/03/21  00:38:05  seidl
 * change sym_libv2.h to chemistry/qc/dmtsym/sym_dmt.h
 *
 * Revision 1.1.1.1  1992/03/17  16:25:54  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:25:53  seidl
 * Initial revision
 *
 * Revision 1.1  1992/02/04  23:48:08  seidl
 * Initial revision
 *
 * Revision 1.6  1992/01/16  19:53:07  seidl
 * use new integral routines
 *
 * Revision 1.5  1992/01/13  19:10:44  seidl
 * remove debug printing, add some comments
 *
 * Revision 1.4  1992/01/09  11:35:59  seidl
 * no longer include util/ipv2/ip_libv2.h
 *
 * Revision 1.3  1991/12/20  16:25:47  seidl
 * add debug print statements, and fix bug in forming scf_vec
 *
 * Revision 1.2  1991/12/17  21:42:21  seidl
 * modified for new util/sgen/sgen
 *
 * Revision 1.1  1991/12/17  14:09:14  seidl
 * Initial revision
 *
 * Revision 1.1  1991/12/02  19:58:51  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#include <stdio.h>
#include <math.h>
#include <tmpl.h>
#include <math/dmt/libdmt.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>

#include "scf.h"

#include "scf_core.gbl"
#include "scf_core.lcl"


GLOBAL_FUNCTION int
scf_core_guess(SCF_VEC,_scf_info,_outfile)
dmt_matrix SCF_VEC;
scf_struct_t *_scf_info;
FILE *_outfile;
{
  int i;
  int errcod;
  int nbasis;
  dmt_matrix S,H,SC,EVECS,SAHALF,SCR;
  double_vector_t evals;

  nbasis = _scf_info->nbfao;
  
  S = dmt_old("libscfv3 overlap matrix");
  H = dmt_old("libscfv3 hcore matrix");

  SC = dmt_columns("libscfv3 scf_core_guess scr1",S);

  EVECS = dmt_create("libscfv3 scf_core_guess scr3",nbasis,COLUMNS);
  SAHALF = dmt_create("libscfv3 scf_core_guess scr4",nbasis,COLUMNS);
  SCR = dmt_create("libscfv3 scf_core_guess scr5",nbasis,COLUMNS);

  errcod = allocbn_double_vector(&evals,"n",nbasis);
  if(errcod != 0) {
    fprintf(_outfile,"scf_core_guess: could not allocate memory for evals\n");
    return(-1);
    }

 /* diagonalize overlap matrix */
 
  dmt_diag(SC,EVECS,evals.d);

 /* form 'sahalf' matrix sahalf = u*ei^-0.5*u~
  * this could be done more efficiently, but what's the point, it takes
  * only 1/10 the time of the two diagonalizations
  */

  for(i=0; i < nbasis ; i++) evals.d[i] = 1.0/sqrt(evals.d[i]);
  dmt_fill(SAHALF,0.0);
  dmt_set_diagonal(SAHALF,evals.d);

  dmt_transpose(EVECS);
  dmt_mult(SAHALF,EVECS,SCR);
  dmt_mult(EVECS,SCR,SAHALF);

 /* now diagonalize core hamiltonian, and multiply vector by SAHALF
  * this way the scf vector will also transform to the S-1/2 basis
  *  c' = S-1/2*c
  *  M(mo) = ~c'*M*c' = ~c*S-1/2*M*S-1/2*c = ~c*M(s-1/2)*c
  */

  dmt_copy(H,SCR);
  dmt_diag(SCR,EVECS,evals.d);
  dmt_mult(SAHALF,EVECS,SCF_VEC);

  free_double_vector(&evals);
  dmt_free(EVECS);
  dmt_free(SCR);
  dmt_free(SC);
  if (!_scf_info->proj_vector) dmt_free(SAHALF);

  return 0;
  }

