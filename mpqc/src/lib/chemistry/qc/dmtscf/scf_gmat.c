
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <tmpl.h>
#include <comm/picl/picl.h>
#include <comm/picl/ext/piclext.h>
#include <math/dmt/libdmt.h>
#include <util/misc/libmisc.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>
#include "scf.h"

#ifndef IOFF
#define IOFF(a,b) ((a)>(b))?(a)*((a)+1)/2+(b):(b)*((b)+1)/2+(a)
#endif

#include "scf_mkj.gbl"
#include "scf_mkk.gbl"
#include "scf_mkgd.gbl"
#include "scf_mkgdlb.gbl"
#include "scf_loopg.gbl"
#include "scf_loopj.gbl"
#include "scf_loopk.gbl"
#include "scf_bnd.gbl"

#include "scf_gmat.gbl"
#include "scf_gmat.lcl"

/* this is a pointer to the buffer which holds the integrals */
double *scf_gmat_intbuf;

/***************************************************************************
 *
 * given a centers struct and an scf struct, initialize some stuff needed
 * to form the G matrix
 *
 */

GLOBAL_FUNCTION int
scf_init_gmat(centers,scf_info)
centers_t *centers;
scf_struct_t *scf_info;
{
  int flags;

  int_initialize_offsets2(centers,centers,centers,centers);

  flags = INT_EREP|INT_NOSTRB|INT_NOSTR1|INT_NOSTR2;

  scf_gmat_intbuf =
    int_initialize_erep(flags,0,centers,centers,centers,centers);
  if (!scf_gmat_intbuf) {
    fprintf(stderr,"scf_init_gmat:  int_initialize_erep() failed\n");
    return -1;
  }

  int_storage(scf_info->int_store);

  if (scf_info->eliminate || scf_info->local_p) {
    if (scf_init_bounds(centers,scf_gmat_intbuf) < 0) {
      fprintf(stderr,"scf_init_gmat:  scf_init_bounds failed\n");
      int_done_erep();
      int_done_offsets2(centers,centers,centers,centers);
      int_done_storage();
      return -1;
    }
  }

  return 0;
}

/*************************************************************************
 *
 * frees memory used by the integral routines and the scf bounds
 */

GLOBAL_FUNCTION VOID
scf_done_gmat(centers, scf_info)
centers_t *centers;
scf_struct_t *scf_info;
{
  int_done_erep();
  int_done_offsets2(centers,centers,centers,centers);
  int_done_storage();

  if (scf_info->eliminate || scf_info->local_p) scf_done_bounds();
}

/**************************************************************************
 *
 * calculates the closed and open shell G matrices
 *
 * input:
 *   scf_info = pointer to initialized scf struct
 *   sym_info = pointer to initialized sym struct
 *   centers  = pointer to initialized centers struct
 *   Gmat     = scattered dmt matrix containing old skeleton G matrix
 *   GmatO    = scattered dmt matrix containing old skeleton G matrix
 *   DPmat    = scattered dmt matrix containing density diff matrix
 *   DPmatO   = scattered dmt matrix containing density diff matrix
 *   SScr1    = scattered dmt scratch matrix
 *   SScr2    = scattered dmt scratch matrix
 *   outfile  = FILE pointer to output (can be null if no output desired)
 *
 * on return:
 *   Gmat and GmatO contain the new skeleton fock matrices
 *
 * return 0 on success and -1 on failure
 */

GLOBAL_FUNCTION int
scf_make_gmat(scf_info,sym_info,centers,Gmat,GmatO,DPmat,DPmatO,
                                                   SScr1,SScr2,outfile)
scf_struct_t *scf_info;
sym_struct_t *sym_info;
centers_t *centers;
dmt_matrix Gmat;
dmt_matrix GmatO;
dmt_matrix DPmat;
dmt_matrix DPmatO;
dmt_matrix SScr1;
dmt_matrix SScr2;
FILE *outfile;
{
  int errcod;

  assert(dmt_distribution(Gmat) == SCATTERED);
  assert(dmt_distribution(DPmat) == SCATTERED);
  assert(dmt_distribution(SScr1) == SCATTERED);
  assert(dmt_distribution(SScr2) == SCATTERED);
  if (scf_info->iopen) {
    assert(dmt_distribution(GmatO) == SCATTERED);
    assert(dmt_distribution(DPmatO) == SCATTERED);
  }

 /* get some timing info */
  tim_enter("gmat");

  if (scf_info->local_p) {
    if (scf_info->load_bal) {
      errcod = scf_make_g_d_lb(centers,scf_info,sym_info,
                               Gmat,GmatO,DPmat,DPmatO,scf_gmat_intbuf,outfile);

    } else if (scf_info->scdft) {
      errcod = scf_make_j_d(centers,scf_info,sym_info,
                               Gmat,DPmat,scf_gmat_intbuf,outfile);

#if 0
     /* this is here until we get the Vxc stuff implemented */
      errcod = scf_make_k_d(centers,scf_info,sym_info,
                               Gmat,DPmat,scf_gmat_intbuf,outfile);
#endif

    } else {
      errcod = scf_make_g_d(centers,scf_info,sym_info,
                               Gmat,GmatO,DPmat,DPmatO,scf_gmat_intbuf,outfile);
    }
  }

 /* if not using local density, then form the g matrix loopwise */
  else {
    if (scf_info->scdft) {
      errcod = scf_make_j_l(centers,scf_info,sym_info,
                            Gmat,DPmat,SScr1,scf_gmat_intbuf,outfile);
      errcod = scf_make_k_l(centers,scf_info,sym_info,
                            Gmat,DPmat,SScr1,scf_gmat_intbuf,outfile);
    } else {
      errcod = scf_make_g_l(centers,scf_info,sym_info,Gmat,GmatO,DPmat,DPmatO,
                            SScr1,SScr2,scf_gmat_intbuf,outfile);
    }
  }

  if (errcod != 0) {
    fprintf(stderr,"scf_gmat: trouble forming gmat 3\n");
    return -1;
  }

  tim_exit("gmat");

  return 0;
}

/***************************************************************************
 *
 * given a dmt density matrix, form an array holding a local copy of the whole
 * matrix.  This is for use by the direct scf functions.
 *
 * input:
 *   scf_info = pointer to initialized scf struct
 *   Pmat     = scattered dmt matrix containing a density matrix 
 *   lp       = pointer to the pointer which will point to the local density
 *              matrix. (uninitialized)
 *   maxp     = pointer to the pointer to the array holding maxP(ishell)
 *   elim     = 1 if doing integral screening
 *
 * on return:
 *   pointer at lp is malloc'd and contains local density
 *   pointer at maxp is malloc'd and contains maxP(ishell)
 *
 * return 0 on success, -1 on failure
 */

GLOBAL_FUNCTION int
scf_make_local_pmat(scf_info,Pmat,lp,maxp,elim)
scf_struct_t *scf_info;
dmt_matrix Pmat;
double **lp;
signed char **maxp;
int elim;
{
  int i,j,ij,ib,jb,isz,jsz,ist,jst,lij;
  int nshell=dmt_nblocks(Pmat);
  int nsht=nshell*(nshell+1)/2;

  double *blk;
  double linv=1.0/log(2.0);
  double tol=pow(2.0,-126.0);
  double tmp,ftmp;
  loop_t *loop;

  assert(dmt_distribution(Pmat) == SCATTERED);

  *lp = (double *) malloc(sizeof(double)*scf_info->nbatri);
  if (!(*lp)) {
    fprintf(stderr,"scf_make_local_pmat:  could not malloc lp\n");
    return -1;
  }

  if (elim) {
    *maxp = (signed char *) malloc(sizeof(signed char)*nsht);
    if (!(*maxp)) {
      fprintf(stderr,"scf_make_local_pmat: could not malloc maxp\n");
      free (*lp);
      return -1;
    }
  }

  loop = dmt_ngl_create("%mr",Pmat);

  while (dmt_ngl_next(loop)) {
    dmt_ngl_create_inner(loop,0);
    while(dmt_ngl_next_inner_m(loop,&ib,&isz,&jb,&jsz,&blk)) {
      dmt_describe_block(Pmat,ib,&ist,&isz);
      dmt_describe_block(Pmat,jb,&jst,&jsz);

      if (ib!=jb) {
        for (i=0; i < isz ; i++) {
          lij=IOFF((i+ist),jst);

          for (j=0; j < jsz ; j++,lij++) {
            (*lp)[lij] = blk[i*jsz+j];
          }
        }
      } else {
        for (i=0; i < isz ; i++) {
          lij=IOFF((i+ist),jst);

          for(j=0; j <= i ; j++,lij++) {
            (*lp)[lij] = blk[i*jsz+j];
          }
        }
      }

      if (elim) {
        tmp=0.0;
        for (i=0; i < isz*jsz; i++) 
          if ((ftmp=fabs(blk[i])) > tmp) tmp=ftmp;

        tmp = (tmp>tol) ? tmp : tol;

        ij = ib*(ib+1)/2+jb;
        (*maxp)[ij] = (signed char) (log(tmp)*linv);
      }
    }
  }

  dmt_ngl_kill(loop);

  return 0;
}

/**************************************************************************
 * 
 * for use by the local_P direct scf routines, this will take an array
 * containing local contributions to the G matrix, and stuff them into a
 * scattered matrix.
 *
 * a global sum should already have been done on lg
 */

GLOBAL_FUNCTION VOID
scf_lgmat_to_scat(lg,Gmat)
double *lg;
dmt_matrix Gmat;
{
  int i,j,ib,jb,isz,jsz,ist,jst,lij;
  int nlocal,nl;
  double *lblk;

  assert(dmt_distribution(Gmat) == SCATTERED);

  nlocal = dmt_nlocal(Gmat);

 /* and transfer to locally held blocks of distributed matrix */
  for (nl=0; nl < nlocal ; nl++) {
    dmt_get_block(Gmat,nl,&ib,&jb,&lblk);
    dmt_describe_block(Gmat,ib,&ist,&isz);
    dmt_describe_block(Gmat,jb,&jst,&jsz);

    if (ib!=jb) {
      for (i=0; i < isz ; i++) {
        lij=IOFF((i+ist),jst);
        for (j=0; j < jsz ; j++,lij++) {
          lblk[i*jsz+j] += lg[lij];
        }
      }
    } else {
      for (i=0; i < isz ; i++) {
        lij=IOFF((i+ist),jst);
        for (j=0; j <= i ; j++,lij++) {
          lblk[i*jsz+j] += lg[lij];
          if(i!=j) lblk[j*jsz+i] += lg[lij];
        }
      }
    }
  }
}

/*****************************************************************************
 *
 * this initialized the arrays needed for the local_P direct fock build.
 *
 * on return:
 *   dpmat, dpmato, maxp, gmat, and gmato are malloc'd and either contain
 *   the local density, or are zero'd
 *
 * return 0 on success, -1 on failure
 */

GLOBAL_FUNCTION int
scf_init_direct_gmat(scf_info, DPmat, DPmatO, dpmat, dpmato, maxp, gmat, gmato)
scf_struct_t *scf_info;
dmt_matrix DPmat;
dmt_matrix DPmatO;
double **dpmat;
double **dpmato;
signed char **maxp;
double **gmat;
double **gmato;
{

  assert(dmt_distribution(DPmat) == SCATTERED);
  if (scf_info->iopen) assert(dmt_distribution(DPmat) == SCATTERED);

  if (scf_make_local_pmat(scf_info,DPmat,dpmat,maxp,scf_info->eliminate) < 0) {
    fprintf(stderr,"scf_init_direct_gmat:  scf_make_local_pmat() failed\n");
    return -1;
  }

  if (scf_info->iopen) {
    if (scf_make_local_pmat(scf_info,DPmatO,dpmato,NULL,0)) {
      fprintf(stderr,"scf_init_direct_gmat:  scf_make_local_pmat() failed\n");

      free(*dpmat);
      if (*maxp) free(*maxp);

      return -1;
    }
  }

 /* now allocate memory for the local G matrices */

  (*gmat) = (double *) malloc(sizeof(double)*scf_info->nbatri);
  if (!(*gmat)) {
    fprintf(stderr,"scf_init_direct_gmat:  could not malloc gmat\n");

    free(*dpmat);
    if (*maxp) free(*maxp);
    if (scf_info->iopen) free(*dpmato);

    return -1;
  }

  bzero(*gmat,sizeof(double)*scf_info->nbatri);

  if (scf_info->iopen) {
    (*gmato) = (double *) malloc(sizeof(double)*scf_info->nbatri);
    if (!(*gmato)) {
      fprintf(stderr,"scf_init_direct_gmat:  could not malloc gmato\n");

      free(*gmat);
      free(*dpmat);
      free(*dpmato);
      if (*maxp) free(*maxp);

      return -1;
    }

  bzero(*gmato,sizeof(double)*scf_info->nbatri);
  }

  return 0;
}

/* free up the arrays used in the local_P direct fock build */

GLOBAL_FUNCTION VOID
scf_done_direct_gmat(dpmat, dpmato, maxp, gmat, gmato)
double *dpmat;
double *dpmato;
signed char *maxp;
double *gmat;
double *gmato;
{
  if (dpmat) free(dpmat);
  if (dpmato) free(dpmato);
  if (maxp) free(maxp);
  if (gmat) free(gmat);
  if (gmato) free(gmato);
}
