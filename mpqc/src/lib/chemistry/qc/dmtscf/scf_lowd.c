
/* $Log$
 * Revision 1.1  1993/12/29 12:53:15  etseidl
 * Initial revision
 *
 * Revision 1.1  1992/07/09  15:44:27  seidl
 * Initial revision
 * */

static char rcsid[]="$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <comm/picl/picl.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <tmpl.h>

#include <math/dmt/libdmt.h>
#include <chemistry/qc/dmtqc/libdmtqc.h>
#include <util/misc/libmisc.h>
#include <util/bio/libbio.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>

#include "scf.h"
#include "scf_lowd.gbl"
#include "scf_lowd.lcl"

/* this function performs a Lowdin analysis.
 * returns: Lowdin bond index in bond_indx
 *          total atomic charge in charge
 */


GLOBAL_FUNCTION int
scf_lowdin(centers,scf_info,irreps,Scf_Vec,bond_indx,charge,outfile)
centers_t *centers;
scf_struct_t *scf_info;
scf_irreps_t *irreps;
dmt_matrix Scf_Vec;
double_matrix_t *bond_indx;
double_vector_t *charge;
FILE *outfile;
{
  int errcod;
  int ndoc=irreps->ir[0].nclosed;
  int nsoc=irreps->ir[0].nopen;
  dmt_matrix S,Pmat,Pmato;

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


 /* form overlap matrix */

  S = dmt_create("overlap",scf_info->nbfao,SCATTERED);

  make_s(S,centers,outfile);


 /* now do the bond index and charges */
  errcod = lowdin_bond_index(Pmat,S,centers,scf_info,bond_indx,charge,outfile);
  if(errcod!=0) {
    fprintf(outfile,"scf_lowdin: lowdin_bond_index failed\n");
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
make_s(S,centers,outfile)
dmt_matrix S;
centers_t *centers;
FILE *outfile;
{
  int i,li,lj,j,isz,ist,jsz,jst;
  int nlocal = dmt_nlocal(S);
  double *data;

  for(i=0; i < nlocal ; i++) {
    dmt_get_block(S,i,&li,&lj,&data);
    int_shell_overlap(centers,centers,data,li,lj);
    }
  }

LOCAL_FUNCTION int
lowdin_bond_index(Pmat,S,centers,scf_info,bindex,qatom,outfile)
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
  dmt_matrix PxS = dmt_create("index mat",centers->nfunc,COLUMNS);
  dmt_matrix PxSt = dmt_create("trans index mat",centers->nfunc,COLUMNS);
  dmt_matrix SC=dmt_columns("column S",S);
  int nlocal=dmt_nlocal(PxS);
  double *pxscol,*pxstcol,*popdata,scale;
  double_vector_t qorb;

  errcod=allocbn_double_vector(&qorb,"n",centers->nfunc);
  if(errcod!=0) return(-1);
  zero_double_vector(&qorb);

  zero_double_vector(qatom);
  zero_double_matrix(bindex);

 /* form S**-1/2 */

  dmt_diag(SC,PxS,qorb.d);

  for(i=0; i < qorb.n; i++) qorb.d[i] = 1.0/sqrt(qorb.d[i]);
  dmt_fill(SC,0.0);
  dmt_set_diagonal(SC,qorb.d);

  dmt_transpose(PxS);
  dmt_mult(SC,PxS,PxSt);
  dmt_mult(PxS,PxSt,SC);

 /* now form S**1/2 = S*S**-1/2 */

  dmt_mult(S,SC,PxSt);

 /* and finally PxS = S**1/2 * P * S**1/2 */
  dmt_mult(Pmat,PxSt,SC);
  dmt_mult(PxSt,SC,PxS);

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

  dmt_copy(PxS,PxSt);
  dmt_transpose(PxSt);

  for(i=0; i < nlocal; i++) {
    dmt_get_col(PxS,i,&j,&pxscol);
    dmt_get_col(PxSt,i,&j,&pxstcol);

    for(j=0; j < qorb.n; j++) pxscol[j] *= pxstcol[j];
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
  dmt_free(SC);
  dmt_free(PxS);

  return(0);
  }
