
#include <stdio.h>
#include <math.h>
#include <tmpl.h>

#include <util/misc/libmisc.h>
#include <util/group/picl.h>
#include <math/dmt/libdmt.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>

#include <chemistry/qc/dmtscf/scf.h>

#include <chemistry/qc/dmtscf/scf_orth.gbl>
#include <chemistry/qc/dmtscf/scf_orth.lcl>

GLOBAL_FUNCTION int
scf_schmidt(scf_info,Scf_Vec,S,all)
scf_struct_t *scf_info;
dmt_matrix Scf_Vec;
dmt_matrix S;
int all;
{
  int ii,i,j;
  int m,ncol;
  int ib,jb,iind,jind,isz,jsz;
  int ihaveit;
  int global_m_ind=0;
  int nbf=scf_info->nbfao;
  int me=mynode0();
  int nlocal=dmt_nlocal(Scf_Vec);
  int nlocalb=dmt_nlocal(S);
  double vtmp;
  double_vector_t cvec,tmpv,scr;
  double *v,*sblk,*ccol;

  tim_enter("schmidt");

  assert(dmt_distribution(Scf_Vec) == COLUMNS);
  assert(dmt_distribution(S) == SCATTERED);

  if (allocbn_double_vector(&tmpv,"n",nbf) != 0) {
    fprintf(stderr,"scf_schmidt: could not alloc tmpv\n");
    return -1;
  }

  if (allocbn_double_vector(&cvec,"n",nbf) != 0) {
    fprintf(stderr,"scf_schmidt: could not alloc cvec\n");
    free_double_vector(&tmpv);
    return(-1);
  }

  if (allocbn_double_vector(&scr,"n",nbf) != 0) {
    fprintf(stderr,"scf_schmidt: could not alloc cvec\n");
    free_double_vector(&tmpv);
    free_double_vector(&cvec);
    return(-1);
  }

  ncol = scf_info->nclosed + scf_info->nopen;
  if(!ncol) ncol++;
  if(all) ncol = nbf;

  while(global_m_ind < ncol) {

  /* find out which node has the column "m" of the vector, and have it
   * broadcast it to all other nodes
   */

    ihaveit=0;
    for(i=0; i < nlocal; i++) {
      dmt_get_col(Scf_Vec,i,&iind,&ccol);
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
      fprintf(stderr,"scf_orth: m = %d ncol = %d\n",m,ncol);
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
      dmt_get_col(Scf_Vec,iind,&i,&ccol);
      if(i>global_m_ind) {
        for(j=0,vtmp=0.0; j<nbf ;j++) vtmp += v[j]*ccol[j];
        for(j=0; j < nbf ; j++) ccol[j] -= vtmp*cvec.d[j];
        }
      }
    
  /* have node responsible for column "m" put it back */

    if(me==ihaveit) {
      for(i=0; i < nlocal; i++) {
        dmt_get_col(Scf_Vec,i,&iind,&ccol);
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
