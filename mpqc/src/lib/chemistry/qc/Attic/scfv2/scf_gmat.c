#define TIME 0

/* $Log$
 * Revision 1.1  1993/12/29 12:53:04  etseidl
 * Initial revision
 *
 * Revision 1.1.1.1  1992/03/17  17:08:43  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:08:41  seidl
 * Initial revision
 *
 * Revision 1.3  1992/01/16  19:53:07  seidl
 * use new integral routines
 *
 * Revision 1.2  1992/01/14  19:33:40  seidl
 * change a few print statements
 *
 * Revision 1.1  1992/01/13  19:13:15  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#include <stdio.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <util/bio/libbio.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/sym/sym_lib.h>
#include "scf.h"

#ifndef IOFF
#define IOFF(a,b) ((a)>(b))?(a)*((a)+1)/2+(b):(b)*((b)+1)/2+(a)
#endif

#include "scf_mkpk.gbl"
#include "scf_mkg.gbl"
#include "scf_mkgd.gbl"

#include "scf_gmat.gbl"
#include "scf_gmat.lcl"

static double *intbuf=NULL;

GLOBAL_FUNCTION int
scf_make_gmat(_scf_info,_sym_info,_irreps,_centers,_mfp,_mast,_cls,_ops,
       _gmat,_gmato,_dpmat,_dpmato,_clbuf,_opbuf,_init,_iter,_outfile)
scf_struct_t *_scf_info;
sym_struct_t *_sym_info;
scf_irreps_t *_irreps;
centers_t *_centers;
scf_mf_ptrs_t *_mfp;
int _mast;
int _cls;
int _ops;
double_vector_t *_gmat;
double_vector_t *_gmato;
sym_d_vector_t *_dpmat;
sym_d_vector_t *_dpmato;
scf_pk_buf_c_t *_clbuf;
scf_pk_buf_o_t *_opbuf;
int _init;
int _iter;
FILE *_outfile;
{
  int i;
  int iopen=_scf_info->iopen;
  int *mptr=&_mfp->cur_ptr;
  int errcod;
  int nbf=_scf_info->nbfao;
  int ntri=_scf_info->nbatri;
  double btime,etime,ttime,clock0();

  double_vector_t ptmp,ptmpo;
  double_matrix_t uaoso,scr;

/* allocate memory for full density matrices, ie not symmetry blocked */

#if TIME
  btime = ttime = clock0();
#endif

  errcod = allocbn_double_vector(&ptmp,"n",ntri);
  if(errcod != 0) {
    fprintf(_outfile,"scf_make_gmat: malloc trouble 1\n");
    return(-1);
    }
  zero_double_vector(&ptmp);
  sym_svec_t_dv(&ptmp,_dpmat);

  if(iopen) {
    errcod = allocbn_double_vector(&ptmpo,"n",ntri);
    if(errcod != 0) {
      fprintf(_outfile,"scf_make_gmat: malloc trouble 2\n");
      return(-1);
      }
    zero_double_vector(&ptmpo);
    sym_svec_t_dv(&ptmpo,_dpmato);
    }

 /* If the density matrices are in the so basis, transform to ao basis using
  * the transform P(ao) = Ua*P(so)*Ua~, where Ua is the aotoso transformation
  * matrix (not the sotoao!).
  *
  * Also note that for this to work, the diagonal elements of P must be
  * multiplied by 2.0 due to the way it is formed in scf_make_density
  */

  if(_scf_info->use_symmetry) {
    init_double_matrix(&uaoso);

    _mfp->cur_ptr = _mfp->uaoso;
    errcod = bread_double_matrix(_mast,&uaoso,&_mfp->cur_ptr);
    if(errcod < 0) {
      fprintf(_outfile,"scf_make_gmat:\n");
      fprintf(_outfile,"trouble reading uaoso from master file\n");
      return(errcod);
      }
    errcod = allocbn_double_matrix(&scr,"n1 n2",nbf,nbf);
    if(errcod != 0) {
      fprintf(_outfile,"scf_make_gmat:\n");
      fprintf(_outfile,"trouble mallocing scr\n");
      return(errcod);
      }

    for(i=0; i < nbf ; i++) ptmp.d[IOFF(i,i)] *= 2.0;
    math_dmxdv_dm(&uaoso,0,&ptmp,0,&scr,0,nbf,nbf,nbf,0);
    math_dmxdm_dv(&scr,0,&uaoso,1,&ptmp,0,nbf,nbf,nbf,0);
    for(i=0; i < nbf ; i++) ptmp.d[IOFF(i,i)] *= 0.5;
    if(iopen) {
      for(i=0; i < nbf ; i++) ptmpo.d[IOFF(i,i)] *= 2.0;
      math_dmxdv_dm(&uaoso,0,&ptmpo,0,&scr,0,nbf,nbf,nbf,0);
      math_dmxdm_dv(&scr,0,&uaoso,1,&ptmpo,0,nbf,nbf,nbf,0);
      for(i=0; i < nbf ; i++) ptmpo.d[IOFF(i,i)] *= 0.5;
      }

    free_double_matrix(&uaoso);
    free_double_matrix(&scr);
    }

#if TIME
  etime = clock0();
  if(mynode0()==0)
    fprintf(_outfile,"scf_gmat: time for pmat = %lf\n",etime-btime);
  btime=etime;
#endif

 /* if this is the first time through, initialize the buffer which will
  * contain the integrals, and then form the needed supermatrix elements
  */
  if(_init) {
    int_initialize_offsets2(_centers,_centers,_centers,_centers);
    intbuf = 
      int_initialize_erep(INT_EREP|INT_NOSTR1|INT_NOSTR2,
       0,_centers,_centers,_centers,_centers);

    int_storage(_scf_info->int_store);

#if TIME
  etime = clock0();
  if(mynode0()==0)
    fprintf(_outfile,"scf_gmat: time to init = %lf\n",etime-btime);
  btime=etime;
#endif

    errcod = scf_make_pk(_centers,_irreps,_scf_info,_sym_info,_mfp,_mast,
           _gmat,_gmato,&ptmp,&ptmpo,_cls,_ops,_clbuf,_opbuf,intbuf,_outfile);
    if(errcod != 0) {
      fprintf(_outfile,"scf_iter: trouble forming pkfiles \n");
      return(-1);
      }
    }
 /* otherwise, just read in supermatrix and form gmat */
  else {
    errcod = scf_make_g(_centers,_irreps,_scf_info,_sym_info,_mfp,_mast,
           _gmat,_gmato,&ptmp,&ptmpo,_cls,_ops,_clbuf,_opbuf,_outfile);
    if(errcod != 0) {
      fprintf(_outfile,"scf_iter: trouble forming gmat 1\n");
      return(-1);
      }
    }

#if TIME
  etime = clock0();
#if 0
  if(mynode0()==0)
#endif
    fprintf(_outfile,"scf_gmat: time for pk = %lf %d\n",etime-btime,mynode0());
  btime=etime;
#endif

 /* finally, calculate any remaining integrals directly and stuff into
  * the g matrix 
  */
  errcod = scf_make_g_d(_centers,_irreps,_scf_info,_sym_info,_mfp,_mast,
             _gmat,_gmato,&ptmp,&ptmpo,intbuf,_iter,_outfile);
  if(errcod != 0) {
    fprintf(_outfile,"scf_iter: trouble forming gmat 2\n");
    return(-1);
    }

  int_reduce_storage_threshold();

#if TIME
  etime = clock0();
#if 0
  if(mynode0()==0)
#endif
    fprintf(_outfile,"scf_gmat: time for gd = %lf  %d\n",etime-btime,mynode0());
  btime=etime;
#endif

  free_double_vector(&ptmp);
  if(iopen) free_double_vector(&ptmpo);

  return 0;
  }
