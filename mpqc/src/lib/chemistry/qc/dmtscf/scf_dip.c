
/* $Log$
 * Revision 1.1  1993/12/29 12:53:15  etseidl
 * Initial revision
 *
 * Revision 1.3  1993/04/27  23:51:45  jannsen
 * fixed a divide by zero problem
 *
 * Revision 1.2  1992/08/12  11:20:01  seidl
 * inlude masses.h from cwd
 *
 * Revision 1.1  1992/07/09  15:44:17  seidl
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
#include "masses.h"

#include "scf.h"
#include "scf_dip.gbl"
#include "scf_dip.lcl"


/* This function calculates and prints the dipole moment wrt the
 * center of mass.  The x, y, and z components of the dipole moment are
 * placed in the vector dipole.
 *
 * It also performs an extended Mulliken analysis (see Huzinaga, et.al.
 * JCP 93(1990)3319),  placing the n*(n+1)/2+n point charges and
 * coordinates for the charges in the expts struct.
 *
 * bond_pops is the output of the function scf_mulliken()
 */

GLOBAL_FUNCTION int
scf_dipole_and_ex_mulliken(centers,scf_info,irreps,Scf_Vec,
        bond_pops,mulpts,dipole,outfile)
centers_t *centers;
scf_struct_t *scf_info;
scf_irreps_t *irreps;
dmt_matrix Scf_Vec;
double_matrix_t *bond_pops;
expts_t *mulpts;
double_vector_t *dipole;
FILE *outfile;
{
  int errcod;
  int ndoc=irreps->ir[0].nclosed;
  int nsoc=irreps->ir[0].nopen;
  dmt_matrix Pmat,Pmato;

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


 /* calculate dipole moment, as well as the charge centers */

  zero_expts(mulpts);

  errcod = dipole_moment(Pmat,centers,scf_info,bond_pops,dipole,mulpts,outfile);
  if(errcod!=0) {
    fprintf(outfile,"scf_dipole_and_ex_mulliken: dipole_moment failed\n");
    return(-1);
    }

 /* free up memory */

  int_done_offsets1(centers,centers);
  int_done_1e();

  dmt_free(Pmat);
  return(0);
  }

LOCAL_FUNCTION int
dipole_moment(Pmat,centers,scf_info,bond_pops,dipole,mulpts,outfile)
dmt_matrix Pmat;
centers_t *centers;
scf_struct_t *scf_info;
double_matrix_t *bond_pops;
double_vector_t *dipole;
expts_t *mulpts;
FILE *outfile;
{
  int i,li,lj,j,isz,ist,jsz,jst;
  int errcod;
  int nat=centers->n;
  int iat,jat,ijat;
  int nfuncmax=int_find_nfuncmax(centers);
  int nlocal=dmt_nlocal(Pmat);

  double dex=0.0,dey=0.0,dez=0.0,dnx=0.0,dny=0.0,dnz=0.0;
  double dtx,dty,dtz,dipmom;
  double *dx=(double *) malloc(sizeof(double)*nfuncmax*nfuncmax);
  double *dy=(double *) malloc(sizeof(double)*nfuncmax*nfuncmax);
  double *dz=(double *) malloc(sizeof(double)*nfuncmax*nfuncmax);
  double *dipder[3];
  double *pmat;
  double scale=1.0;
  double debye=2.54176548;
  double com[3];

  com[0] = com[1] = com[2] = 0.0;

  if (nfuncmax!=0 && (dx==NULL || dy == NULL || dz == NULL)) {
    fprintf(outfile,"scf_dipole_moment: malloc failed\n");
    return(-1);
    }

  dipder[0] = dx;
  dipder[1] = dy;
  dipder[2] = dz;

 /* calculate the center of mass */

  errcod=center_of_mass(centers,scf_info,com,outfile);
  if(errcod!=0) return(-1);
  

 /* form dipole contributions 
  *
  *  dex = sum_i sum_j Pij*<i|x|j>
  *
  * also form contributions to the center of the overlap population
  *
  * <x>ab = sum_i(in a) sum_j(in b) Pij*<i|x|j>  
  *         -----------------------------------
  *           sum_i(in a) sum_j(in b) Pij*Sij
  */

  for(i=0; i < nlocal ; i++) {
    dmt_get_block(Pmat,i,&li,&lj,&pmat);

    int_shell_dipole(centers,centers,com,dipder,li,lj);

    dmt_describe_block(Pmat,li,&ist,&isz);
    dmt_describe_block(Pmat,lj,&jst,&jsz);

    scale=(li==lj)?1.0:2.0;

    for(j=0; j < isz*jsz; j++) {
      dex+=dx[j]*pmat[j]*scale;
      dey+=dy[j]*pmat[j]*scale;
      dez+=dz[j]*pmat[j]*scale;
      }

    iat=centers->center_num[li];
    jat=centers->center_num[lj];

    scale=(iat==jat && li!=lj) ? 2.0 : 1.0;
    for(j=0; j < isz*jsz; j++) {
      ijat=(iat>jat)? iat*(iat+1)/2+jat : jat*(jat+1)/2+iat;
      mulpts->p[ijat].r[0] += scale*pmat[j]*dx[j];
      mulpts->p[ijat].r[1] += scale*pmat[j]*dy[j];
      mulpts->p[ijat].r[2] += scale*pmat[j]*dz[j];
      }
    }

  /* sum the contributions to the center of overlap population */

  for(iat=ijat=0; iat < nat; iat++) {
    for(jat=0; jat <= iat; jat++,ijat++) {
      scale=(iat==jat) ? 1.0 : 2.0;
      gsum0(mulpts->p[ijat].r,3,5,mtype_get(),0);
      bcast0(mulpts->p[ijat].r,sizeof(double)*3,mtype_get(),0);

      if (fabs(bond_pops->d[iat][jat]) >= 1.0e-15) {
         mulpts->p[ijat].r[0] /= bond_pops->d[iat][jat];
         mulpts->p[ijat].r[1] /= bond_pops->d[iat][jat];
         mulpts->p[ijat].r[2] /= bond_pops->d[iat][jat];
         }
      mulpts->p[ijat].charge = -bond_pops->d[iat][jat]*scale;
      }
    }

  /* sum dipole contributions */

  gsum0(&dex,1,5,mtype_get(),0);
  gsum0(&dey,1,5,mtype_get(),0);
  gsum0(&dez,1,5,mtype_get(),0);

 /* form nuclear contributions to dipole moment */

  ijat=nat*(nat+1)/2;
  for(i=0; i < nat; i++) {
    dnx += centers->center[i].charge*(centers->center[i].r[0]-com[0]);
    dny += centers->center[i].charge*(centers->center[i].r[1]-com[1]);
    dnz += centers->center[i].charge*(centers->center[i].r[2]-com[2]);

    mulpts->p[i+ijat].charge=centers->center[i].charge;
    mulpts->p[i+ijat].r[0]=centers->center[i].r[0]-com[0];
    mulpts->p[i+ijat].r[1]=centers->center[i].r[1]-com[1];
    mulpts->p[i+ijat].r[2]=centers->center[i].r[2]-com[2];
    }

  dipole->d[0] = dnx-dex;
  dipole->d[1] = dny-dey;
  dipole->d[2] = dnz-dez;

  dex*= -debye;
  dey*= -debye;
  dez*= -debye;
  dnx*= debye;
  dny*= debye;
  dnz*= debye;

  dtx=dex+dnx;
  dty=dey+dny;
  dtz=dez+dnz;

  if(mynode0()==0) {
    fprintf(outfile,"\n  dipole moments wrt center of mass in debye\n");
    fprintf(outfile,"\n  dex = %12.8f dey = %12.8f dez = %12.8f\n",dex,dey,dez);
    fprintf(outfile,"  dnx = %12.8f dny = %12.8f dnz = %12.8f\n",dnx,dny,dnz);
    fprintf(outfile,"  dtx = %12.8f dty = %12.8f dtz = %12.8f\n",dtx,dty,dtz);

    dipmom = sqrt(dtx*dtx+dty*dty+dtz*dtz);
    fprintf(outfile,
      "\n  total dipole moment wrt center of mass = %15.10f debye\n",dipmom);
    }

  free(dx);
  free(dy);
  free(dz);

  return(0);
  }

LOCAL_FUNCTION int
center_of_mass(centers,scf_info,com,outfile)
centers_t *centers;
scf_struct_t *scf_info;
double *com;
FILE *outfile;
{
  int i;
  double X,Y,Z,M;
  double *x=(double *) malloc(sizeof(double)*centers->n);
  double *y=(double *) malloc(sizeof(double)*centers->n);
  double *z=(double *) malloc(sizeof(double)*centers->n);
  double *m=(double *) malloc(sizeof(double)*centers->n);

  if (centers->n!=0 && (x==NULL || y == NULL || z == NULL || m == NULL)) {
    fprintf(outfile,"center_of_mass: malloc failed\n");
    return(-1);
    }

  X=Y=Z=M=0.0;
  for(i=0; i < centers->n; i++) {
    x[i]=centers->center[i].r[0];
    y[i]=centers->center[i].r[1];
    z[i]=centers->center[i].r[2];
    m[i]=Atomic_Mass[(int)centers->center[i].charge];
    M+=m[i];
    X+=m[i]*x[i];
    Y+=m[i]*y[i];
    Z+=m[i]*z[i];
    }

  X /= M;
  Y /= M;
  Z /= M;

  com[0]=X;
  com[1]=Y;
  com[2]=Z;

  if(mynode0()==0) {
    int j; double num;
    double xx,yy,zz;
    double ausq_to_angsq = 0.2800283608302436;
    double A, B, C;
    double h = 6.626176e-34;
    double c = 2.99792458e10;
    double pisq = 9.869604404;
    double conv = 1.660565e-47;
    double cm_to_hz = 2.9979243e4;

    double_vector_t inert,evals;
    double_matrix_t evecs;

    allocbn_double_vector(&inert,"n",6);
    zero_double_vector(&inert);
    allocbn_double_vector(&evals,"n",3);
    allocbn_double_matrix(&evecs,"n1 n2",3,3);

    fprintf(outfile,"\n  center of mass in a.u.\n    %f %f %f\n\n",X,Y,Z);
  
    for (i=0; i < centers->n ; i++) {
      x[i] -= X;
      y[i] -= Y;
      z[i] -= Z;

      inert.d[0] += ausq_to_angsq*m[i]*(y[i]*y[i] + z[i]*z[i]);
      inert.d[1] -= ausq_to_angsq*m[i]*x[i]*y[i];
      inert.d[2] += ausq_to_angsq*m[i]*(x[i]*x[i] + z[i]*z[i]);
      inert.d[3] -= ausq_to_angsq*m[i]*x[i]*z[i];
      inert.d[4] -= ausq_to_angsq*m[i]*y[i]*z[i];
      inert.d[5] += ausq_to_angsq*m[i]*(x[i]*x[i] + y[i]*y[i]);
      }

    math_diag_dv(&inert,&evals,&evecs,1,1.0e-15,1);

    fprintf(outfile,"\n  Principal moments of inertia (amu*A**2)\n");
    fprintf(outfile,"       Ia = %19.10f\n",evals.d[0]);
    fprintf(outfile,"       Ib = %19.10f\n",evals.d[1]);
    fprintf(outfile,"       Ic = %19.10f\n",evals.d[2]);

    A = h/(8.0*pisq*evals.d[0]*c*conv);
    B = h/(8.0*pisq*evals.d[1]*c*conv);
    C = h/(8.0*pisq*evals.d[2]*c*conv);

    fprintf(outfile,"\n  Rotational constants\n");
    fprintf(outfile,"\n             cm**-1         MHz\n");
    fprintf(outfile,"  A = %19.10f %19.10f\n",A,A*cm_to_hz);
    fprintf(outfile,"  B = %19.10f %19.10f\n",B,B*cm_to_hz);
    fprintf(outfile,"  C = %19.10f %19.10f\n",C,C*cm_to_hz);

    for (i=0; i < centers->n ; i++) {
      for (j=0; j < 3; j++) {
        num = 0.0;
        num += evecs.d[0][j]*x[i];
        num += evecs.d[1][j]*y[i];
        num += evecs.d[2][j]*z[i];
        switch (j) {
          case 0:
            xx = num;
            break;
          case 1:
            yy = num;
            break;
          case 2:
            zz = num;
          }
        }
      x[i]=xx;
      y[i]=yy;
      z[i]=zz;
      }

    free_double_vector(&inert);
    free_double_vector(&evals);
    free_double_matrix(&evecs);
    }

  free(x);
  free(y);
  free(z);
  free(m);

  return(0);
  }
