
/* $Log$
 * Revision 1.1  1993/12/29 12:53:16  etseidl
 * Initial revision
 *
 * Revision 1.1  1992/07/09  15:44:31  seidl
 * Initial revision
 * */

static char rcsid[]="$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <tmpl.h>

#include <comm/picl/picl.h>

#include <math/dmt/libdmt.h>

#include <util/misc/libmisc.h>
#include <util/bio/libbio.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/dmtqc/libdmtqc.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>

#include "scf.h"
#include "scf_mull.gbl"
#include "scf_mull.lcl"

/* this function performs a Mulliken analysis.
 * returns: overlap bond populations in bond_pops
 *          Mulliken bond index in bond_indx
 *          total atomic charge in charge
 */

GLOBAL_FUNCTION int
scf_mulliken(centers,scf_info,irreps,Scf_Vec,bond_pops,bond_indx,charge,outfile)
centers_t *centers;
scf_struct_t *scf_info;
scf_irreps_t *irreps;
dmt_matrix Scf_Vec;
double_matrix_t *bond_pops;
double_matrix_t *bond_indx;
double_vector_t *charge;
FILE *outfile;
{
  int errcod;
  int ndoc=irreps->ir[0].nclosed;
  int nsoc=irreps->ir[0].nopen;
  dmt_matrix S,Pmat,Pmato,Pop;

 /* initialize centers struct */

  int_initialize_1e(0,0,centers,centers);
  int_initialize_offsets1(centers,centers);

 /* form density matrices */

  Pmat = dmt_create("density",scf_info->nbfao,SCATTERED);
  if(scf_info->iopen)
    Pmato = dmt_create("open density",scf_info->nbfao,SCATTERED);

  dmt_density(Scf_Vec,ndoc,Pmat);
  if(scf_info->iopen) {
    dmt_open_density(Scf_Vec,ndoc,nsoc,Pmato);
    dmt_sum_scaled(Pmato,0.5,Pmat);
    dmt_free(Pmato);
    }
  dmt_scale(Pmat,2.0);


 /* form overlap matrix and AO populations */

  S = dmt_create("overlap",scf_info->nbfao,SCATTERED);
  Pop = dmt_create("population",scf_info->nbfao,SCATTERED);

  make_s_and_pop(S,Pmat,Pop,centers,outfile);


 /* calculate overlap bond populations */

  zero_double_matrix(bond_pops);

  errcod=bond_pop(Pop,centers,bond_pops,outfile);
  if(errcod!=0) {
    fprintf(outfile,"scf_mulliken: bond_pop failed\n");
    return(-1);
    }

  dmt_free(Pop);


 /* now do the bond index and charges */

  errcod=mulliken_bond_index(Pmat,S,centers,scf_info,bond_indx,charge,outfile);
  if(errcod!=0) {
    fprintf(outfile,"scf_mulliken: mulliken_bond_index failed\n");
    return(-1);
    }

 /* free up memory */

  int_done_offsets1(centers,centers);
  int_done_1e();

  dmt_free(S);
  dmt_free(Pmat);
  return(0);
  }


LOCAL_FUNCTION VOID
make_s_and_pop(S,Pmat,Pop,centers,outfile)
dmt_matrix S;
dmt_matrix Pmat;
dmt_matrix Pop;
centers_t *centers;
FILE *outfile;
{
  int i,li,lj,j,isz,ist,jsz,jst;
  int nlocal = dmt_nlocal(S);
  double *data,*pdata,*popdata;

 /* Pop_ij = S_ij * P_ij */

  for(i=0; i < nlocal ; i++) {
    dmt_get_block(S,i,&li,&lj,&data);
    int_shell_overlap(centers,centers,data,li,lj);

    dmt_get_block(Pmat,i,&li,&lj,&pdata);
    dmt_get_block(Pop,i,&li,&lj,&popdata);
    dmt_describe_block(Pmat,li,&ist,&isz);
    dmt_describe_block(Pmat,lj,&jst,&jsz);

    for(j=0; j < isz*jsz; j++) popdata[j] = data[j]*pdata[j];
    }
  }

LOCAL_FUNCTION int
bond_pop(Pop,centers,batom,outfile)
dmt_matrix Pop;
centers_t *centers;
double_matrix_t *batom;
FILE *outfile;
{
  int i,li,lj,j,isz,ist,jsz,jst;
  int iat,jat;
  int nat=centers->n;
  int nlocal = dmt_nlocal(Pop);
  double *popdata,scale;


 /* the overlap population Q is
  *
  * Qab = sum_i(in a) sum_j (in b) Sij*Pij
  */

  for(i=0; i < nlocal ; i++) {
    dmt_get_block(Pop,i,&li,&lj,&popdata);
    dmt_describe_block(Pop,li,&ist,&isz);
    dmt_describe_block(Pop,lj,&jst,&jsz);

    iat=centers->center_num[li];
    jat=centers->center_num[lj];

    scale=(iat==jat && li!=lj) ? 2.0 : 1.0;
    for(j=0; j < isz*jsz; j++) batom->d[iat][jat] += scale*popdata[j];

    batom->d[jat][iat]=batom->d[iat][jat];
    }

  popdata=(double *) malloc(sizeof(double)*batom->n2);
  if(popdata==NULL) return(-1);

  for(i=0; i < batom->n1; i++) {
    gop1(batom->d[i],batom->n2,popdata,'+',mtype_get());
    }

  free(popdata);

  return(0);
  }

LOCAL_FUNCTION int
mulliken_bond_index(Pmat,S,centers,scf_info,bindex,qatom,outfile)
dmt_matrix Pmat;
dmt_matrix S;
centers_t *centers;
scf_struct_t *scf_info;
double_matrix_t *bindex;
double_vector_t *qatom;
FILE *outfile;
{
  int i,j,errcod;
  int li,lj,isz,ist,jsz,jst,iat,jat;
  double *pxscol,*pxstcol,*popdata,scale;
  dmt_matrix PxS = dmt_create("index mat",centers->nfunc,COLUMNS);
  dmt_matrix PxSt = dmt_create("trans index mat",centers->nfunc,COLUMNS);
  double_vector_t qorb;
  int nlocal=dmt_nlocal(PxS);


 /* form matrix PxS = PS */

  dmt_mult(Pmat,S,PxS);

  dmt_copy(PxS,PxSt);
  dmt_transpose(PxSt);


 /* initialize arrays */

  errcod=allocbn_double_vector(&qorb,"n",centers->nfunc);
  if(errcod!=0) return(-1);

  zero_double_vector(&qorb);

  zero_double_vector(qatom);
  zero_double_matrix(bindex);


 /* the trace of PxS == n_elec
  * the charge on atom a, Qa = sum_i(in a) PxSaa
  */

  dmt_get_diagonal(PxS,qorb.d);

  for(i=0; i < centers->nshell; i++) {
    int at=centers->center_num[i];
    int last=(i==centers->nshell-1)?centers->nfunc:centers->func_num[i+1];
    for(j=centers->func_num[i]; j < last; j++) qatom->d[at] -= qorb.d[j];
    }
      
  for(i=0; i < centers->n; i++) qatom->d[i] += centers->center[i].charge;


 /* now for the bond index.  don't ask me, I just copied it from
  * bondex
  *
  * bond index ab = sum_i(in a) sum_j(in b) PxSij * PxSji
  */

  for(i=0; i < nlocal; i++) {
    dmt_get_col(PxS,i,&j,&pxscol);
    dmt_get_col(PxSt,i,&j,&pxstcol);

    for(j=0; j < centers->nfunc; j++) pxscol[j] *= pxstcol[j];
    }
    
  dmt_free(PxSt);
  PxSt = dmt_scatter("scat index",PxS);

  nlocal=dmt_nlocal(PxSt);

  for(i=0; i < nlocal ; i++) {
    dmt_get_block(PxSt,i,&li,&lj,&popdata);
    dmt_describe_block(PxSt,li,&ist,&isz);
    dmt_describe_block(PxSt,lj,&jst,&jsz);

    iat=centers->center_num[li];
    jat=centers->center_num[lj];

    scale = (iat==jat && li!=lj) ? 2.0 : 1.0;

    for(j=0; j < isz*jsz; j++) bindex->d[iat][jat] += scale*popdata[j];

    bindex->d[jat][iat]=bindex->d[iat][jat];
    }

  for(i=0; i < bindex->n1; i++) {
    gop1(bindex->d[i],bindex->n2,qorb.d,'+',mtype_get());
    }

  free_double_vector(&qorb);
  dmt_free(PxS);
  return(0);
  }
