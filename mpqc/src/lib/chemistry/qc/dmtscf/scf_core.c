
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <tmpl.h>

#include <util/misc/libmisc.h>
#include <math/dmt/libdmt.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>


#include <chemistry/qc/dmtscf/scf.h>

#include <chemistry/qc/dmtscf/scf_core.gbl>
#include <chemistry/qc/dmtscf/scf_core.lcl>

/***********************************************************************
 * 
 * given a dmt matrix and an scf struct, construct a trial scf vector 
 * using the core hamiltonian guess.
 *
 * input:
 *   scf_info = pointer to initialized scf struct
 *   Scf_Vec  = column distributed matrix
 *   Hcore    = scattered dmt matrix containing AO core hamiltonian
 *   S        = scattered dmt matrix containing AO overlap integrals
 *   Sahalf   = column distributed matrix
 *
 * on return:
 *   Scf_Vec points to the trial vector
 *   Sahalf contains the S**-1/2 matrix
 *
 * return 0 on success, -1 on failure
 */

GLOBAL_FUNCTION int
scf_core_guess(scf_info,Scf_Vec, Hcore, S, Sahalf)
scf_struct_t *scf_info;
dmt_matrix Scf_Vec;
dmt_matrix Hcore;
dmt_matrix S;
dmt_matrix Sahalf;
{
  int i;
  int nbasis;
  dmt_matrix SC,Evecs,Scr;
  double *evals;

  assert(dmt_distribution(Scf_Vec) == COLUMNS);
  assert(dmt_distribution(Hcore) == SCATTERED);
  assert(dmt_distribution(S) == SCATTERED);
  assert(dmt_distribution(Sahalf) == COLUMNS);

  nbasis = scf_info->nbfao;
  
  SC = dmt_columns("dmtscf scf_core_guess scr1",S);

  Evecs = dmt_create("dmtscf scf_core_guess scr3",nbasis,COLUMNS);
  Scr = dmt_create("dmtscf scf_core_guess scr5",nbasis,COLUMNS);

  evals = (double *) malloc(sizeof(double)*nbasis);
  if (!evals) {
    fprintf(stderr,"scf_core_guess:  could not malloc evals(%d)\n",nbasis);
    return -1;
  }

 /* diagonalize overlap matrix */
 
  dmt_diag(SC,Evecs,evals);

 /* form 'sahalf' matrix sahalf = u*ei^-0.5*u~
  * this could be done more efficiently, but what's the point, it takes
  * only 1/10 the time of the two diagonalizations
  */

  for (i=0; i < nbasis ; i++) evals[i] = 1.0/sqrt(evals[i]);
  dmt_fill(Sahalf,0.0);
  dmt_set_diagonal(Sahalf,evals);

  dmt_transpose(Evecs);
  dmt_mult(Sahalf,Evecs,Scr);
  dmt_mult(Evecs,Scr,Sahalf);

 /* now diagonalize core hamiltonian, and multiply vector by SAHALF
  * this way the scf vector will also transform to the S-1/2 basis
  *  c' = S-1/2*c
  *  M(mo) = ~c'*M*c' = ~c*S-1/2*M*S-1/2*c = ~c*M(s-1/2)*c
  */

#if 0
  dmt_mult(Hcore,Sahalf,Evecs); /* Evecs is a scratch matrix here. */
  dmt_mult(Sahalf,Evecs,Scr);   /* Evecs is a scratch matrix here. */
  dmt_diag(Scr,Evecs,evals);
#else
  {
    dmt_matrix Hc = dmt_columns("foo",Hcore);
    dmt_diag(Hc,Evecs,evals);
    dmt_free(Hc);
    if (scf_info->print_flg&2048) {
      if (mynode0()==0) printf("HCore eigenvectors");
      dmt_printf("%15.7f ",Evecs);
      if (mynode0()==0) {
        printf("eigenvalues\n");
        for (i=0; i < nbasis; i++)
          printf("%5d %20.10f\n",i+1,evals[i]);
      }
    }
  }
#endif
  dmt_mult(Sahalf,Evecs,Scf_Vec);

  free(evals);

  dmt_free(Evecs);
  dmt_free(Scr);
  dmt_free(SC);

  return 0;
}
