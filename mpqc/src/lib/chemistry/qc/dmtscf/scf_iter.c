
/* $Log$
 * Revision 1.1  1993/12/29 12:53:15  etseidl
 * Initial revision
 *
 * Revision 1.16  1992/06/29  17:48:08  seidl
 * use scf_schmidt() for non-local_p as well
 *
 * Revision 1.15  1992/06/23  20:04:32  seidl
 * change dmt matrices to uppercase,
 * get rid of unnecessary matrice
 *
 * Revision 1.14  1992/06/17  21:54:18  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.13  1992/06/16  16:27:08  seidl
 * reset first three iterations for cheat
 *
 * Revision 1.12  1992/05/19  20:58:30  seidl
 * add cheat stuff
 *
 * Revision 1.11  1992/05/04  11:06:38  seidl
 * remove pk ints, timings controlled by print_flg now
 *
 * Revision 1.10  1992/04/17  14:56:57  seidl
 * pass overlap matrix to scf_schmit()
 *
 * Revision 1.9  1992/04/13  11:06:56  seidl
 * exit after printing out energy if maxiter exceeded
 *
 * Revision 1.8  1992/04/08  20:39:43  seidl
 * move tim_print statement so times are printed out after first iteration
 *
 * Revision 1.7  1992/04/06  19:34:11  seidl
 * use iter%ckpt_freq, not iter%5
 *
 * Revision 1.6  1992/04/06  12:36:36  seidl
 * reset converged flag and diis_er
 *
 * Revision 1.5  1992/04/01  01:03:37  seidl
 * fix bounds
 *
 * Revision 1.4  1992/03/31  22:25:43  seidl
 * print timing info if _scf_info->print_tim is true
 *
 * Revision 1.3  1992/03/21  00:39:17  seidl
 * change sym_libv2.h to chemistry/qc/dmtsym/sym_dmt.h
 *
 * Revision 1.2  1992/03/18  12:16:14  seidl
 * add timing statement
 *
 * Revision 1.1.1.1  1992/03/17  16:26:23  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:26:22  seidl
 * Initial revision
 *
 * Revision 1.7  1992/02/26  12:21:46  seidl
 * free fdiag when done with it
 *
 * Revision 1.6  1992/02/21  15:17:58  seidl
 * add schmit orthogonalization
 *
 * Revision 1.5  1992/02/13  00:42:31  seidl
 * free up memory in diis when done
 * get correct form of mo fock matrices for open shell when done
 *
 * Revision 1.4  1992/02/10  17:01:29  seidl
 * remove clock0() call
 *
 * Revision 1.3  1992/02/07  12:59:59  seidl
 * remove timing stuff, don't checkpoint after 1st iteration
 *
 * Revision 1.2  1992/02/05  14:05:23  seidl
 * checkpoint vector and fock matrices every 5 iterations
 * use util/misc/libmisc timing routines
 *
 * Revision 1.1  1992/02/04  23:48:08  seidl
 * Initial revision
 *
 * Revision 1.8  1992/01/16  19:53:07  seidl
 * use new integral routines
 *
 * Revision 1.7  1992/01/13  19:14:26  seidl
 * add some timing statements, call scf_make_gmat instead of
 * scf_make_g_d and scf_make_g
 *
 * Revision 1.6  1992/01/09  11:58:50  seidl
 * first parallel version
 *
 * Revision 1.5  1992/01/06  11:48:04  seidl
 * add call to scf_schmit
 *
 * Revision 1.4  1992/01/02  16:20:13  seidl
 * add open-shell and diis
 *
 * Revision 1.3  1991/12/24  19:30:41  seidl
 * fix bug in writing of open-shell gmatrix to scr array
 *
 * Revision 1.2  1991/12/24  11:47:47  seidl
 * clean up matrices when converged
 *
 * Revision 1.1  1991/12/20  16:23:08  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#include <stdio.h>
#include <tmpl.h>
#include <comm/picl/picl.h>
#include <math/dmt/libdmt.h>
#include <util/misc/libmisc.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>
#include "scf.h"

#include "scf_mkden.gbl"
#include "scf_gmat.gbl"
#include "scf_en.gbl"
#include "scf_fock.gbl"
#include "scf_diis.gbl"
#include "scf_orth.gbl"

#include "scf_iter.gbl"
#include "scf_iter.lcl"

GLOBAL_FUNCTION int
scf_iter(_scf_info,_sym_info,_irreps,_centers,
       SCF_VEC,FOCK,FOCKO,_evals,_occ_num,
       GMAT,GMATO,PMAT,DPMAT,PMATO,DPMATO,SCR1,SCR2,SCR3,SSCR1,SSCR2,_outfile)
scf_struct_t *_scf_info;
sym_struct_t *_sym_info;
scf_irreps_t *_irreps;
centers_t *_centers;
dmt_matrix SCF_VEC;
dmt_matrix FOCK;
dmt_matrix FOCKO;
double_vector_t *_evals;
double_vector_t *_occ_num;
dmt_matrix GMAT;
dmt_matrix GMATO;
dmt_matrix PMAT;
dmt_matrix DPMAT;
dmt_matrix PMATO;
dmt_matrix DPMATO;
dmt_matrix SCR1;
dmt_matrix SCR2;
dmt_matrix SCR3;
dmt_matrix SSCR1;
dmt_matrix SSCR2;
FILE *_outfile;
{
  int i;
  int iter,maxiter;
  int iopen=_scf_info->iopen;
  int errcod;
  int nbasis=_scf_info->nbfao;
  char vecfile[512],fockfile[512],fockofile[512];

  dmt_matrix S=dmt_old("libscfv3 overlap matrix");

  maxiter = _scf_info->maxiter;

/* what are the names of the checkpoint files? */

  sprintf(vecfile,"%s%s.scfvec",_scf_info->ckptdir,_scf_info->fname);
  sprintf(fockfile,"%s%s.fock",_scf_info->ckptdir,_scf_info->fname);
  sprintf(fockofile,"%s%s.focko",_scf_info->ckptdir,_scf_info->fname);

/* begin the iteration here */

  _scf_info->converged=0;
  _scf_info->diis_er=0.0;
  for(iter=0; iter < maxiter ; iter++) {

   /* beware, there are side effects.  on exit SSCR2 contains the full
    * open-shell part of the fock matrix 
    */
    make_fock(_centers,_scf_info,_sym_info,FOCK,FOCKO,GMAT,GMATO,SSCR1,SSCR2,_outfile);

    scf_calculate_energy(_scf_info,FOCK,PMAT,DPMAT,
                                   SSCR2,PMATO,SSCR1,iter,_outfile);

  /* perform diis extrapolation, returns new fock matrix in ao basis */
    if(_scf_info->diis_flg) {
      errcod = scf_diis(_scf_info,FOCK,FOCKO,SCF_VEC,_occ_num,iter,
              SCR1,SCR2,SCR3,0,_outfile);
      if(errcod != 0) {
        fprintf(_outfile,"scf_iter: trouble in scf_diis\n");
        return(-1);
        }
      }

  /* if open-shell, form effective fock matrix used to get new vector */
    if(iopen) {
      errcod = scf_make_fock(_scf_info,FOCK,FOCKO,SCF_VEC,_occ_num,
                    SCR1,SCR2,SCR3,iter,_outfile);
      if(errcod != 0) {
        fprintf(_outfile,"scf_iter: trouble in scf_make_fock\n");
        return(-1);
        }
      }

/* transform fock matrix to mo basis, then diagonalize to obtain new vector
 * scf_make_fock already does this, so don't do it for open shell
 */

    if(!iopen) {
      dmt_mult(FOCK,SCF_VEC,SCR1);
      dmt_mult(SCF_VEC,SCR1,SCR2);
      dmt_copy(SCR2,FOCK);
      }

    dmt_copy(FOCK,SCR2); /* dmt_diag needs a columns dist. matrix */
    dmt_diag(SCR2,SCR1,_evals->d);
    dmt_transpose(SCF_VEC);
    dmt_mult(SCF_VEC,SCR1,SCR2);
    dmt_copy(SCR2,SCF_VEC);

/* un-level shift eigenvalues */

    if(iopen) {
      double_vector_t fdiag;
      double occi;
      double occ0=_occ_num->d[0];

      errcod = allocbn_double_vector(&fdiag,"n",nbasis);
      dmt_get_diagonal(FOCK,fdiag.d);
       
      for(i=0; i < nbasis; i++) {
        occi=_occ_num->d[i];
        if(occi==occ0 && occi) {
          _evals->d[i] += _scf_info->lvl_shift;
          fdiag.d[i] += _scf_info->lvl_shift;
          }
        else if(occi) {
          _evals->d[i] += 0.5*_scf_info->lvl_shift;
          fdiag.d[i] += 0.5*_scf_info->lvl_shift;
          }
        }
      dmt_set_diagonal(FOCK,fdiag.d);
    
      free_double_vector(&fdiag);
      }

    if(_scf_info->print_flg & 2) {
      sync0();
      tim_print(0);
      }

/* if converged, return */

    if(_scf_info->converged) {
 /* make sure the real closed and open shell mo fock matrices are placed in
  * FOCK and FOCKO.
  */
      if(_scf_info->diis_flg) {
        scf_diis(_scf_info,FOCK,FOCKO,SCF_VEC,_occ_num,iter,
              SCR1,SCR2,SCR3,1,_outfile);
        }
      else {
        make_fock(_centers,_scf_info,_sym_info,
                  FOCK,FOCKO,GMAT,GMATO,SSCR1,SSCR2,_outfile);

        dmt_mult(FOCK,SCF_VEC,SCR1);
        dmt_mult(SCF_VEC,SCR1,SCR2);
        dmt_copy(SCR2,FOCK);
        if(iopen) {
          dmt_mult(FOCKO,SCF_VEC,SCR1);
          dmt_mult(SCF_VEC,SCR1,SCR2);
          dmt_copy(SCR2,FOCKO);
          }
        }

      if(mynode0()==0) fprintf(_outfile,"\n  converged\n");
      return(0);
      }

/* orthogonalize new vector */

    errcod = scf_schmidt(_scf_info,_irreps,SCF_VEC,S,1,_outfile);
    if(errcod!=0) {
      fprintf(_outfile,"trouble orthogonalizing scf vector\n");
      return(-1);
      }

/* make new density matrices */

    errcod = scf_make_density(_scf_info,_irreps,
      SCF_VEC,PMAT,DPMAT,PMATO,DPMATO,_occ_num,_outfile);
    if(errcod!=0) {
      fprintf(_outfile,"trouble forming density matrices\n");
      return(-1);
      }

/* reset density if appropriate */

    if(_scf_info->eliminate && 
      ((iter+1)%_scf_info->p_reset_freq == 0 || 
       (_scf_info->cheat && iter < 3))) {
      dmt_copy(PMAT,DPMAT);
      dmt_fill(GMAT,0.0);
      if(_scf_info->iopen) {
        dmt_copy(PMATO,DPMATO);
        dmt_fill(GMATO,0.0);
        }
      if(mynode0()==0) fprintf(_outfile,"  resetting density matrices\n");
      }

/* and armed with the new density matrix, form new gmatrix */

    if((iter+1) < maxiter) {
      errcod = scf_make_gmat(_scf_info,_sym_info,_irreps,_centers,
               GMAT,GMATO,DPMAT,DPMATO,SSCR1,SSCR2,0,iter+1,_outfile);
      if(errcod != 0) {
        fprintf(_outfile,"scf_iter: trouble forming gmat\n");
        return(-1);
        }
      }

/* checkpoint every "ckpt_freq" iterations */
    if((iter+1)%_scf_info->ckpt_freq==0) {
      if(mynode0()==0)
        fprintf(_outfile,"  checkpointing vector and fock matrices\n");
      dmt_write(vecfile,SCF_VEC);
      dmt_write(fockfile,FOCK);
      if(iopen) dmt_write(fockofile,FOCKO);
      }

    } /* ad nauseum */


/* clean up after yourself */
  if(_scf_info->diis_flg) {
    scf_diis(_scf_info,FOCK,FOCKO,SCF_VEC,_occ_num,iter,
          SCR1,SCR2,SCR3,1,_outfile);
    }
  else {
    make_fock(_centers,_scf_info,_sym_info,FOCK,FOCKO,GMAT,GMATO,SSCR1,SSCR2,_outfile);

    dmt_mult(FOCK,SCF_VEC,SCR1);
    dmt_mult(SCF_VEC,SCR1,SCR2);
    dmt_copy(SCR2,FOCK);
    if(iopen) {
      dmt_mult(FOCKO,SCF_VEC,SCR1);
      dmt_mult(SCF_VEC,SCR1,SCR2);
      dmt_copy(SCR2,FOCKO);
      }
    }

  return 0;
  }

LOCAL_FUNCTION VOID
make_fock(_centers,_scf_info,_sym_info,FOCK,FOCKO,GMAT,GMATO,SSCR1,SSCR2,_outfile)
centers_t *_centers;
scf_struct_t *_scf_info;
sym_struct_t *_sym_info;
dmt_matrix FOCK;
dmt_matrix FOCKO;
dmt_matrix GMAT;
dmt_matrix GMATO;
dmt_matrix SSCR1;
dmt_matrix SSCR2;
FILE *_outfile;
{
  int iopen=_scf_info->iopen;
  dmt_matrix HCORE=dmt_old("libscfv3 hcore matrix");

 /* form full g matrix from the skeleton gmatrix, place the result in SSCR1 */

  if(_sym_info->g > 1) {
    sym_sym_matrix(_centers,_sym_info,GMAT,SSCR1,_outfile);
    if(iopen) sym_sym_matrix(_centers,_sym_info,GMATO,SSCR2,_outfile);
    }
  else {
    dmt_copy(GMAT,SSCR1);
    if(iopen) dmt_copy(GMATO,SSCR2);
    }

 /* F = H + G
  * FO = H + G - GO */

  dmt_copy(HCORE,FOCK);
  dmt_sum(SSCR1,FOCK);

  if(iopen) {
    dmt_copy(SSCR2,FOCKO);
    dmt_scale(FOCKO,-1.0);
    dmt_sum(FOCK,FOCKO);
    }
  }
