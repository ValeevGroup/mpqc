
/* $Log$
 * Revision 1.1  1993/12/29 12:53:14  etseidl
 * Initial revision
 *
 * Revision 1.15  1992/06/30  11:57:03  seidl
 * don't do orthogonalization again if projected vector is used
 *
 * Revision 1.14  1992/06/29  17:48:59  seidl
 * use scf_schmidt() for non-local_p as well
 *
 * Revision 1.13  1992/06/23  20:04:45  seidl
 * change dmt matrices to uppercase,
 * get rid of unnecessary matrice
 *
 * Revision 1.12  1992/06/17  21:54:33  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.11  1992/06/16  16:24:22  seidl
 * add tim_print() call
 *
 * Revision 1.10  1992/05/13  18:23:24  jannsen
 * Added fflush calls.
 *
 * Revision 1.9  1992/05/04  11:09:14  seidl
 * remove pk integral stuff
 *
 * Revision 1.8  1992/04/22  16:01:14  seidl
 * add oldvecfile for vector projection stuff, print out eigenvalues if
 * debug = 1
 *
 * Revision 1.7  1992/04/17  14:58:00  seidl
 * pass overlap matrix to scf_schmit
 *
 * Revision 1.6  1992/04/07  18:09:15  jannsen
 * doesn't write the final matrices, lets mpqcnode do this so the user
 * can control whether or not they are written
 *
 * Revision 1.5  1992/04/07  18:04:20  jannsen
 *
 * Revision 1.4  1992/04/06  12:38:55  seidl
 * time oeis, don't read anything in if restart
 *
 * Revision 1.3  1992/04/01  01:04:45  seidl
 * fix bounds checking
 *
 * Revision 1.2  1992/03/21  00:41:23  seidl
 * change sym_libv2.h to chemistry/qc/dmtsym/sym_dmt.h
 *
 * Revision 1.1.1.1  1992/03/17  16:26:55  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:26:54  seidl
 * Initial revision
 *
 * Revision 1.10  1992/03/04  15:59:47  seidl
 * call scf_done_bounds
 *
 * Revision 1.9  1992/02/26  17:19:15  seidl
 * free some buffers for scf_mkg
 *
 * Revision 1.8  1992/02/26  12:40:51  seidl
 * have the user create vector and fock matrices, no longer pass in pointers
 * to them
 * free oei matrices
 *
 * Revision 1.7  1992/02/26  12:22:28  seidl
 * free pk buffers when done
 *
 * Revision 1.6  1992/02/21  15:17:43  seidl
 * add schmit orthogonalization
 *
 * Revision 1.5  1992/02/18  16:17:25  seidl
 * allocate memory for integral matrices here
 *
 * Revision 1.4  1992/02/13  00:43:27  seidl
 * don't write matrices to local directory, let the
 * user do that instead
 *
 * Revision 1.3  1992/02/10  17:07:25  seidl
 * remove timing calls
 *
 * Revision 1.2  1992/02/05  13:59:35  seidl
 * add restart and checkpoint stuff, use util/misc/libmisc timing routines
 *
 * Revision 1.1  1992/02/04  23:48:08  seidl
 * Initial revision
 *
 * Revision 1.10  1992/01/16  19:49:40  seidl
 * use new integral routines, call int_done_x when cleaning up
 *
 * Revision 1.9  1992/01/13  19:17:07  seidl
 * add timing statements, replace scf_make_g_d and scf_make_pk with
 * scf_make_gmat
 *
 * Revision 1.8  1992/01/09  14:42:25  seidl
 * print statement causes icc to choke
 *
 * Revision 1.7  1992/01/09  11:53:10  seidl
 * add parallel code
 *
 * Revision 1.6  1992/01/02  18:12:30  seidl
 * have scf_file place total scf energy from the last calculation in
 * the current scf_struct
 *
 * Revision 1.5  1992/01/02  16:23:28  seidl
 * a few cosmetic changes
 *
 * Revision 1.4  1991/12/24  19:34:03  seidl
 * use new zero_ routines, get rid of zerov and zerom
 *
 * Revision 1.3  1991/12/24  11:49:41  seidl
 * add pretty print stuff
 *
 * Revision 1.2  1991/12/20  16:32:19  seidl
 * finish initial revision
 *
 * Revision 1.1  1991/12/17  21:45:25  seidl
 * Initial revision
 *
 * Revision 1.1  1991/12/02  19:58:51  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#include <stdio.h>
#include <tmpl.h>
#include <math/dmt/libdmt.h>
#include <util/misc/libmisc.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>

#include "scf.h"

#include "scf_oeis.gbl"
#include "scf_core.gbl"
#include "scf_mkden.gbl"
#include "scf_gmat.gbl"
#include "scf_iter.gbl"
#include "scf_orth.gbl"
#include "scf_bnd.gbl"
#include "scf_ex.gbl"

#include "scf_vect.gbl"
#include "scf_vect.lcl"

GLOBAL_FUNCTION int
scf_vector(_scf_info,_sym_info,_irreps,_centers,
           FOCK,FOCKO,SCF_VEC,_outfile)
scf_struct_t *_scf_info;
sym_struct_t *_sym_info;
scf_irreps_t *_irreps;
centers_t *_centers;
dmt_matrix FOCK;
dmt_matrix FOCKO;
dmt_matrix SCF_VEC;
FILE *_outfile;
{
  int j;
  int errcod;
  int nbasis=_scf_info->nbfao;
  char vecfile[512],fockfile[512],fockofile[512],oldvecfile[512];

  dmt_matrix PMAT,PMATO,DPMAT,DPMATO,GMAT,GMATO;
  dmt_matrix S,T,V,H;
  dmt_matrix SCR1,SCR2,SCR3,SSCR1,SSCR2;
  double_vector_t occ_num,evals;

/* what are the names of the checkpoint files? */

  sprintf(vecfile,"%s%s.scfvec",_scf_info->ckptdir,_scf_info->fname);
  sprintf(oldvecfile,"%s%s.oldvec",_scf_info->ckptdir,_scf_info->fname);
  sprintf(fockfile,"%s%s.fock",_scf_info->ckptdir,_scf_info->fname);
  sprintf(fockofile,"%s%s.focko",_scf_info->ckptdir,_scf_info->fname);

/* calculate one-electron integrals and place them in the master file */

  S = dmt_create("libscfv3 overlap matrix",nbasis,SCATTERED);
  T = dmt_create("libscfv3 kinetic matrix",nbasis,SCATTERED);
  V = dmt_create("libscfv3 potential matrix",nbasis,SCATTERED);
  H = dmt_create("libscfv3 hcore matrix",nbasis,SCATTERED);

  tim_enter("scf_oeis");
  errcod = scf_oeis(_scf_info,_sym_info,_irreps,_centers,S,T,V,H,_outfile);
  tim_exit("scf_oeis");
  if(errcod != 0) {
    fprintf(_outfile,"scf_vector:\n");
    fprintf(_outfile,"trouble in forming one-electron integrals\n");
    return(-1);
    }

/* we don't need the T and V matrices any longer, so free up the memory */
  dmt_free(T);
  dmt_free(V);

/* get old vector if there is one, otherwise, construct a guess */

  if(_scf_info->restart) {
    if(mynode0()==0) fprintf(_outfile,"\n  using old vector\n\n");
    }
  else if(_scf_info->warmrestart) {
    dmt_read(vecfile,SCF_VEC);
    if(mynode0()==0) fprintf(_outfile,
      "\n  read vector from checkpoint file %s\n\n",vecfile);
    }
  else if(_scf_info->proj_vector) {
    tim_enter("proj_vector");
    if(mynode0()==0) fprintf(_outfile,
      "\n  forming projection of old scf vector\n\n");
    errcod = scf_core_guess(SCF_VEC,_scf_info,_outfile);
    if(errcod != 0) {
      fprintf(_outfile,"scf_vector:\n");
      fprintf(_outfile,"trouble forming guess scf vector\n");
      return(-1);
      }
    errcod = scf_project_vector(_centers,_irreps,_scf_info,SCF_VEC,
      oldvecfile,_outfile);
    tim_exit("proj_vector");
    if(errcod != 0) {
      fprintf(_outfile,"scf_vector:\n");
      fprintf(_outfile,"trouble forming projected guess scf vector\n");
      return(-1);
      }
    }
  else {
    if(_scf_info->print_flg & 16) tim_print(0);
    if(mynode0()==0) fprintf(_outfile,
      "\n  first run, so defaulting to core-hamiltonian guess\n\n");
    errcod = scf_core_guess(SCF_VEC,_scf_info,_outfile);
    if(errcod != 0) {
      fprintf(_outfile,"scf_vector:\n");
      fprintf(_outfile,"trouble forming guess scf vector\n");
      return(-1);
      }
    }
  if(mynode0()==0) fflush(_outfile);

/* set up occupation numbers and initialize eigenvalues */

  errcod = allocbn_double_vector(&evals,"n",nbasis);
  if(errcod != 0) {
    fprintf(_outfile,"could not allocate memory for evals vector\n");
    return(-1);
    }
  errcod = allocbn_double_vector(&occ_num,"n",nbasis);
  if(errcod != 0) {
    fprintf(_outfile,"could not allocate memory for occ_num vector\n");
    return(-1);
    }

  for(j=0; j < _irreps->ir[0].nclosed ; j++)                occ_num.d[j]=2.0;
  for(; j < _irreps->ir[0].nclosed+_irreps->ir[0].nopen ; j++) 
                                                            occ_num.d[j]=1.0;
  for(; j < _irreps->ir[0].num_so ; j++)                    occ_num.d[j]=0.0;

/* form initial density matrices */

  PMAT = dmt_create("libscfv3 density matrix",nbasis,SCATTERED);
  dmt_fill(PMAT,0.0);
  DPMAT = dmt_create("libscfv3 density diff matrix",nbasis,SCATTERED);
  dmt_fill(DPMAT,0.0);
  if(_scf_info->iopen) {
    PMATO = dmt_create("libscfv3 open density matrix",nbasis,SCATTERED);
    dmt_fill(PMATO,0.0);
    DPMATO = dmt_create("libscfv3 open density diff matrix",nbasis,SCATTERED);
    dmt_fill(DPMATO,0.0);
    }
  else {
    PMATO = dmt_nil();
    DPMATO = dmt_nil();
    }

/* orthogonalize vector unless a projected guess was used, in which case
 * the orthogonalization has already been done
 */

  if(!_scf_info->proj_vector || _scf_info->restart || _scf_info->warmrestart) {
    errcod = scf_schmidt(_scf_info,_irreps,SCF_VEC,S,0,_outfile);
    if(errcod!=0) {
      fprintf(_outfile,"trouble orthogonalizing vector\n");
      return(-1);
      }
    }

  errcod = scf_make_density(_scf_info,_irreps,
                           SCF_VEC,PMAT,DPMAT,PMATO,DPMATO,&occ_num,_outfile);
  if(errcod!=0) {
    fprintf(_outfile,"trouble forming density matrices\n");
    return(-1);
    }

/* and let's form the pk file */

  GMAT = dmt_create("libscfv3 g matrix",nbasis,SCATTERED);
  dmt_fill(GMAT,0.0);
  if(_scf_info->iopen) {
    GMATO = dmt_create("libscfv3 open g matrix",nbasis,SCATTERED);
    dmt_fill(GMATO,0.0);
    }
  else {
    GMATO = dmt_nil();
    }

/* allocate scratch arrays now */

  SCR1 = dmt_create("scf_iter: scr1",nbasis,COLUMNS);
  SCR2 = dmt_create("scf_iter: scr2",nbasis,COLUMNS);
  SCR3 = dmt_create("scf_iter: scr3",nbasis,COLUMNS);
  SSCR1 = dmt_create("scf_iter: scr4",nbasis,SCATTERED);
  SSCR2 = dmt_create("scf_iter: scr5",nbasis,SCATTERED);

/* form g matrix */
  errcod = scf_make_gmat(_scf_info,_sym_info,_irreps,_centers,
                         GMAT,GMATO,DPMAT,DPMATO,SSCR1,SSCR2,1,0,_outfile);
  if(errcod != 0) {
    fprintf(_outfile,"scf_vector: trouble forming gmat\n");
    return(-1);
    }

/* now iterate */

  errcod = scf_iter(_scf_info,_sym_info,_irreps,_centers,
             SCF_VEC,FOCK,FOCKO,&evals,&occ_num,
             GMAT,GMATO,PMAT,DPMAT,PMATO,DPMATO,SCR1,SCR2,SCR3,SSCR1,SSCR2,
             _outfile);
  if(errcod != 0) {
    fprintf(_outfile,"scf_vector: trouble in scf_iter\n");
    return(-1);
    }

/* release memory used for two-electron integral routines */
  int_done_erep();
  int_done_offsets2(_centers,_centers,_centers,_centers);
  int_done_storage();
  if(_scf_info->eliminate || _scf_info->local_p) scf_done_bounds();

/* compute the exchange energy, if requested */
  if (_scf_info->exchange) {
      scf_ex(_scf_info,_centers,PMAT);
    }

#if FINAL_CHECKPOINT /* This can be done by mpqcnode. */
/* final checkpoint, just in case... */

  if(mynode0()==0) {
    fprintf(_outfile,"  final checkpoint\n");
    fflush(_outfile);
    }
  dmt_write(vecfile,SCF_VEC);
  dmt_write(fockfile,FOCK);
  if(_scf_info->iopen) dmt_write(fockofile,FOCKO);
#endif

/* delete checkpoint files if desired */
  if(_scf_info->ckpt_del && (mynode0()==0)) {
    sprintf(vecfile,"%s%s.scfvec",_scf_info->ckptdir,_scf_info->fname);
    sprintf(fockfile,"%s%s.fock",_scf_info->ckptdir,_scf_info->fname);
    sprintf(fockofile,"%s%s.focko",_scf_info->ckptdir,_scf_info->fname);

    unlink(vecfile);
    unlink(fockfile);
    unlink(fockofile);
    if(mynode0()==0) fprintf(_outfile,"  deleted checkpoint files\n");
    }

/* print eigenvalues if desired */

  if((_scf_info->print_flg&1) && mynode0()==0) {
    print_evals(_outfile,&evals,&occ_num);
    }

/* clean up your room young man */
  dmt_free(S);
  dmt_free(H);
  dmt_free(PMAT);
  dmt_free(GMAT);
  dmt_free(DPMAT);
  dmt_free(SCR1);
  dmt_free(SCR2);
  dmt_free(SCR3);
  dmt_free(SSCR1);
  dmt_free(SSCR2);
  if(_scf_info->iopen) {
    dmt_free(PMATO);
    dmt_free(GMATO);
    dmt_free(DPMATO);
    }
  free_double_vector(&occ_num);
  free_double_vector(&evals);

  return(0);
  }

LOCAL_FUNCTION VOID
print_evals(_outfile,evals,occ_num)
FILE *_outfile;
double_vector_t *evals;
double_vector_t *occ_num;
{
  int i,j;

  fprintf(_outfile,"\n  eigenvalues and occupation numbers\n");

  for (i=i=0; i < evals->n; i+=6) {
    for(j=0; (j<6) && (i+j<evals->n); j++) fprintf(_outfile,"%9d    ",i+j);
    fprintf(_outfile,"\n");
    for(j=0; (j<6) && (i+j<evals->n); j++) 
      fprintf(_outfile,"%13.6f",evals->d[i+j]);
    fprintf(_outfile,"\n");
    for(j=0; (j<6) && (i+j<evals->n); j++) 
      fprintf(_outfile,"%13.6f",occ_num->d[i+j]);
    fprintf(_outfile,"\n\n");
    }
  }
