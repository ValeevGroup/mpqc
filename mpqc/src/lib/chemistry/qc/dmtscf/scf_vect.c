
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <tmpl.h>
#include <comm/picl/picl.h>
#include <comm/picl/ext/piclext.h>
#include <math/dmt/libdmt.h>
#include <util/misc/libmisc.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>

#include "scf.h"

#include "scf_oeis.gbl"
#include "scf_core.gbl"
#include "scf_iter.gbl"
#include "scf_orth.gbl"
#include "scf_proj.gbl"
#include "scf_mkden.gbl"
#include "scf_ex.gbl"

#include "scf_vect.gbl"
#include "scf_vect.lcl"

/************************************************************************
 *
 * this function will return a converged scf vector
 *
 * input:
 *   scf_info   = pointer to scf struct
 *   sym_info   = pointer to sym struct
 *   centers    = pointer to centers struct
 *   Fock       = scattered dmt matrix (overwritten)
 *   FockO      = scattered dmt matrix (overwritten)
 *   Scf_Vec    = columns dmt matrix (overwritten)
 *   oldcenters = pointer to center struct using old basis (can be null)
 *                only used for vector projection
 *   outfile    = FILE pointer to output
 *
 * on return:
 *   Scf_Vec contains the converged scf_vector
 *   Fock and FockO contain the AO basis fock matrices
 *   scf_info->converged is set to 1 if the calculation did indeed converge
 *
 * return 0 on success, -1 on failure
 */

GLOBAL_FUNCTION int
scf_vector(scf_info,sym_info,centers,Fock,FockO,Scf_Vec,oldcenters,outfile)
scf_struct_t *scf_info;
sym_struct_t *sym_info;
centers_t *centers;
dmt_matrix Fock;
dmt_matrix FockO;
dmt_matrix Scf_Vec;
centers_t *oldcenters;
FILE *outfile;
{
  int j;
  int nbasis=scf_info->nbfao;
  char vecfile[512],fockfile[512],fockofile[512],oldvecfile[512];

  dmt_matrix S,T,V,Hcore;
  double *occ_num, *evals;

  assert(dmt_distribution(Scf_Vec) == COLUMNS);
  assert(dmt_distribution(Fock) == SCATTERED);
  if (scf_info->iopen) assert(dmt_distribution(FockO) == SCATTERED);

 /* what are the names of the checkpoint files? */

  sprintf(vecfile,"%s%s.scfvec",scf_info->ckptdir,scf_info->fname);
  sprintf(fockfile,"%s%s.fock",scf_info->ckptdir,scf_info->fname);
  sprintf(fockofile,"%s%s.focko",scf_info->ckptdir,scf_info->fname);

 /* the oldvec file should be in the current directory */

  sprintf(oldvecfile,"./%s.oldvec",scf_info->fname);

 /* calculate one-electron integrals */

  S = dmt_create("dmtscf overlap matrix",nbasis,SCATTERED);
  T = dmt_create("dmtscf kinetic matrix",nbasis,SCATTERED);
  V = dmt_create("dmtscf potential matrix",nbasis,SCATTERED);
  Hcore = dmt_create("dmtscf hcore matrix",nbasis,SCATTERED);

  if (scf_oeis(scf_info,sym_info,centers,S,T,V,Hcore,outfile) < 0) {
    fprintf(stderr,"scf_vector:  trouble forming one-electron integrals\n");
    return -1;
  }

 /* we don't need the T and V matrices any longer, so free up the memory */
  dmt_free(T);
  dmt_free(V);

 /* get old vector if there is one, otherwise, construct a guess */

  if (scf_info->restart) {
    if (mynode0()==0 && outfile)
      fprintf(outfile,"\n  scf_vect: using old vector\n\n");

  } else if (scf_info->warmrestart) {
    dmt_read(vecfile,Scf_Vec);
    if (mynode0()==0 && outfile) {
      fprintf(outfile,"\n  scf_vect: read vector from checkpoint file %s\n\n",
                      vecfile);
    }

  } else if (scf_info->proj_vector) {
    dmt_matrix Sahalf;

    tim_enter("proj_vector");
    if (mynode0()==0 && outfile)
      fprintf(outfile,"\n  scf_vect: forming projection of old scf vector\n\n");

    Sahalf = dmt_create("scf_vect scratch matrix",nbasis,COLUMNS);

    if (scf_core_guess(scf_info,Scf_Vec,Hcore,S,Sahalf) < 0) {
      fprintf(stderr,"scf_vector:  trouble forming guess scf vector\n");
      return -1;
    }

    if (scf_project_vector(centers,scf_info,Scf_Vec,S,Sahalf,oldvecfile,
                                oldcenters,outfile) < 0) {
      fprintf(stderr,"scf_vector: "
                     "trouble forming projected guess scf vector\n");
      return -1;
    }

    dmt_free(Sahalf);
    tim_exit("proj_vector");

  } else {
    dmt_matrix Sahalf;

    if (scf_info->print_flg & 16) tim_print(0);

    if (mynode0()==0 && outfile) {
      fprintf(outfile,"\n");
      fprintf(outfile,"  first run, so defaulting to core-hamiltonian guess");
      fprintf(outfile,"\n\n");
    }

    Sahalf = dmt_create("scf_vect scratch matrix",nbasis,COLUMNS);

    if (scf_core_guess(scf_info,Scf_Vec,Hcore,S,Sahalf) < 0) {
      fprintf(stderr,"scf_vector: trouble forming guess scf vector\n");
      return -1;
    }

    dmt_free(Sahalf);
  }

  if (mynode0()==0 && outfile) fflush(outfile);

/* set up occupation numbers and initialize eigenvalues */

  evals = (double *) malloc(sizeof(double)*nbasis);
  if (!evals) {
    fprintf(stderr,"scf_vect: could not allocate memory for evals\n");
    return -1;
  }

  occ_num = (double *) malloc(sizeof(double)*nbasis);
  if (!occ_num) {
    fprintf(stderr,"scf_vect: could not allocate memory for occ_num\n");
    return -1;
  }

  for (j=0; j < scf_info->nclosed ; j++)                occ_num[j]=2.0;
  for (; j < scf_info->nclosed+scf_info->nopen ; j++)   occ_num[j]=1.0;
  for (; j < nbasis ; j++)                              occ_num[j]=0.0;

 /* orthogonalize vector unless a projected guess was used, in which case
  * the orthogonalization has already been done
  */

  if (!scf_info->proj_vector || scf_info->restart || scf_info->warmrestart) {
    if (scf_schmidt(scf_info,Scf_Vec,S,0) < 0) {
      fprintf(stderr,"scf_vector:  trouble orthogonalizing vector\n");
      return -1;
    }
  }

 /* now iterate */

  if (scf_do_iter(scf_info,sym_info,centers,Scf_Vec,Hcore,S,Fock,FockO,
                       evals,occ_num,outfile) < 0) {
    fprintf(stderr,"scf_vector: trouble in scf_iter\n");
    return -1;
  }

 /* compute the exchange energy, if requested */
  if (scf_info->exchange) {
    scf_ex(scf_info,centers,Scf_Vec);
    if (mynode0()==0 && outfile) {
      fprintf(outfile,"\n  The exchange energy is %14.10f\n",scf_info->e_exc);
    }
  }

 /* delete checkpoint files if desired */
  if (scf_info->ckpt_del && (mynode0()==0)) {
    sprintf(vecfile,"%s%s.scfvec",scf_info->ckptdir,scf_info->fname);
    sprintf(fockfile,"%s%s.fock",scf_info->ckptdir,scf_info->fname);
    sprintf(fockofile,"%s%s.focko",scf_info->ckptdir,scf_info->fname);

    unlink(vecfile);
    unlink(fockfile);
    unlink(fockofile);
    if (mynode0()==0 && outfile)
      fprintf(outfile,"  deleted checkpoint files\n");
  }

 /* print eigenvalues if desired */

  if ((scf_info->print_flg&1) && mynode0()==0 && outfile) {
    print_evals(outfile,evals,occ_num,nbasis);
  }

/* clean up your room young man */
  dmt_free(S);
  dmt_free(Hcore);

  free(occ_num);
  free(evals);

  return 0;
}

LOCAL_FUNCTION VOID
print_evals(outfile,evals,occ_num,n)
FILE *outfile;
double *evals;
double *occ_num;
int n;
{
  int i,j;

  fprintf(outfile,"\n  eigenvalues and occupation numbers\n");

  for (i=i=0; i < n; i+=6) {
    for (j=0; (j<6) && (i+j<n); j++) fprintf(outfile,"%9d    ",i+j);
    fprintf(outfile,"\n");

    for (j=0; (j<6) && (i+j<n); j++) 
      fprintf(outfile,"%13.6f",evals[i+j]);
    fprintf(outfile,"\n");

    for (j=0; (j<6) && (i+j<n); j++) 
      fprintf(outfile,"%13.6f",occ_num[i+j]);
    fprintf(outfile,"\n\n");
  }
}
