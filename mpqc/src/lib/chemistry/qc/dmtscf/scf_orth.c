/* $Log$
 * Revision 1.1  1993/12/29 12:53:16  etseidl
 * Initial revision
 *
 * Revision 1.10  1993/04/27  23:54:25  jannsen
 * changed mynode() to mynode0()
 *
 * Revision 1.9  1992/06/30  11:52:31  seidl
 * free cmat before receiving it
 *
 * Revision 1.8  1992/06/29  17:49:37  seidl
 * use dmt matrices now
 *
 * Revision 1.7  1992/06/23  20:04:42  seidl
 * change dmt matrices to uppercase,
 * get rid of unnecessary matrice
 *
 * Revision 1.6  1992/06/17  21:53:10  jannsen
 * cleaned up for saber-c and changed to ngl loops
 *
 * Revision 1.5  1992/05/26  20:18:09  jannsen
 * use mtype_get to get message types for global operations
 * check results of memory allocations
 *
 * Revision 1.4  1992/04/22  16:00:03  seidl
 * fix bug in dm_to_col
 *
 * Revision 1.3  1992/04/17  14:57:37  seidl
 * overlap matrix passed in now
 *
 * Revision 1.2  1992/03/21  00:41:09  seidl
 * change sym_libv2.h to chemistry/qc/dmtsym/sym_dmt.h
 *
 * Revision 1.1.1.1  1992/03/17  16:26:53  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:26:52  seidl
 * Initial revision
 *
 * Revision 1.2  1992/02/21  15:16:35  seidl
 * modify for dmt, only use if local_p
 *
 * Revision 1.1  1992/02/04  23:48:08  seidl
 * Initial revision
 *
 * Revision 1.2  1992/01/16  19:53:07  seidl
 * use new integral routines
 *
 * Revision 1.1  1992/01/03  12:44:49  seidl
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

#include "scf_orth.gbl"
#include "scf_orth.lcl"

GLOBAL_FUNCTION int
scf_schmidt(_scf_info,_irreps,SCF_VEC,S,all,_outfile)
scf_struct_t *_scf_info;
scf_irreps_t *_irreps;
dmt_matrix SCF_VEC;
dmt_matrix S;
int all;
FILE *_outfile;
{
  int ii,i,j,jlst;
  int m,ncol;
  int errcod;
  int ib,jb,iind,jind,isz,jsz;
  int ihaveit;
  int global_m_ind=0;
  int nbf=_scf_info->nbfao;
  int me=mynode0();
  int dim=cubedim0();
  int nlocal=dmt_nlocal(SCF_VEC);
  int nlocalb=dmt_nlocal(S);
  double vtmp,scale;
  double_vector_t cvec,tmpv,scr;
  double *v,*sblk,*ccol;

  tim_enter("schmidt");

  errcod = allocbn_double_vector(&tmpv,"n",nbf);
  if(errcod != 0) {
    fprintf(_outfile,"scf_schmidt: could not alloc tmpv\n");
    return(-1);
    }

  errcod = allocbn_double_vector(&cvec,"n",nbf);
  if(errcod != 0) {
    fprintf(_outfile,"scf_schmidt: could not alloc cvec\n");
    free_double_vector(&tmpv);
    return(-1);
    }

  errcod = allocbn_double_vector(&scr,"n",nbf);
  if(errcod != 0) {
    fprintf(_outfile,"scf_schmidt: could not alloc cvec\n");
    free_double_vector(&tmpv);
    free_double_vector(&cvec);
    return(-1);
    }

  ncol = _irreps->ir[0].nclosed + _irreps->ir[0].nopen;
  if(!ncol) ncol++;
  if(all) ncol = nbf;

  while(global_m_ind < ncol) {

  /* find out which node has the column "m" of the vector, and have it
   * broadcast it to all other nodes
   */

    ihaveit=0;
    for(i=0; i < nlocal; i++) {
      dmt_get_col(SCF_VEC,i,&iind,&ccol);
      if(iind==global_m_ind) {
        ihaveit=me;
        for(j=0; j < cvec.n; j++) cvec.d[j]=ccol[j];
        }
      }
    gsum0(&ihaveit,1,2,mtype_get(),0);
    bcast0(&ihaveit,sizeof(int),mtype_get(),0);
   
  /* bcast0_double_vector will do a malloc every time, so free cvec */

    if(me!=ihaveit) free_double_vector(&cvec);
    bcast0_double_vector(&cvec,0,ihaveit);
      
  /* now that each node has a copy of the current column, do the nasty */
    zero_double_vector(&tmpv);
    v=tmpv.d;

  /* v[i] = sum_j C[j][m]*S[i][j] */

    for(ii=0; ii < nlocalb; ii++) {
      dmt_get_block(S,ii,&ib,&jb,&sblk);
      dmt_describe_block(S,ib,&iind,&isz);
      dmt_describe_block(S,jb,&jind,&jsz);
      if(ib==jb) {
        for(i=0; i < isz ; i++) {
          for(j=0; j < jsz ; j++) {
            v[i+iind] += cvec.d[j+jind]*sblk[i*jsz+j];
            }
          }
        }
      else {
        for(i=0; i < isz ; i++) {
          for(j=0; j < jsz ; j++) {
            v[i+iind] += cvec.d[j+jind]*sblk[i*jsz+j];
            v[j+jind] += cvec.d[i+iind]*sblk[i*jsz+j];
            }
          }
        }
      }

    gop1(tmpv.d,tmpv.n,scr.d,'+',mtype_get());
    v=tmpv.d;

  /* vtmp = sum_i v[i]*c[i][m] 
   *      = sum_i sum_j (c[j][m]*s[i][j])*c[i][m]
   */

    for(i=0,vtmp=0.0; i < nbf ; i++) vtmp += v[i]*cvec.d[i];
    if(vtmp == 0.0) {
      fprintf(_outfile,"scf_orth: m = %d ncol = %d\n",m,ncol);
      continue;
      }

    if(vtmp < 1.0e-20) vtmp = 1.0e-20;
    vtmp = 1.0/sqrt(vtmp);

  /* v[i] = v[i]/sqrt(vtmp)
   * c[i][m] = c[i][m]/sqrt(vtmp)
   */

    for(i=0; i < nbf ; i++) {
      v[i] *= vtmp;
      cvec.d[i] *= vtmp;
      }
  /*                             _
   * vtmp = sum_j v[j]*c[j][i]    |
   * c[j][i] -= vtmp*c[j][m]     _| for all (i>m)
   */

    for(iind=0; iind<nlocal; iind++) {
      dmt_get_col(SCF_VEC,iind,&i,&ccol);
      if(i>global_m_ind) {
        for(j=0,vtmp=0.0; j<nbf ;j++) vtmp += v[j]*ccol[j];
        for(j=0; j < nbf ; j++) ccol[j] -= vtmp*cvec.d[j];
        }
      }
    
  /* have node responsible for column "m" put it back */

    if(me==ihaveit) {
      for(i=0; i < nlocal; i++) {
        dmt_get_col(SCF_VEC,i,&iind,&ccol);
        if(iind==global_m_ind) {
          for(j=0; j < cvec.n; j++) ccol[j]=cvec.d[j];
          }
        }
      }

    global_m_ind++;
    }

  free_double_vector(&cvec);
  free_double_vector(&tmpv);
  free_double_vector(&scr);
  
  tim_exit("schmidt");

  return(0);
  }
