
/* $Log$
 * Revision 1.1  1993/12/29 12:53:15  etseidl
 * Initial revision
 *
 * Revision 1.2  1992/08/12  11:20:31  seidl
 * include vdw_rad.h from cwd
 *
 * Revision 1.1  1992/07/09  15:44:25  seidl
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
#include "vdw_rad.h"


#include "scf.h"
#include "scf_esp.gbl"
#include "scf_esp.lcl"

#ifndef M_PI
#define M_PI    3.14159265358979323846
#define M_PI_2  1.57079632679489661923
#endif

/* This function is given the results of an extended mulliken analysis
 * in mulpts, and the x, y, and z components of the dipole moment (in au)
 * in dipole.  The points in mulpts are used to calculate a classical
 * electrostatic potential, which is then used to calculate atomic
 * point charges.  The charges are returned in charges.
 *
 * See Tasi, Kiricsi, and Forster, J.Comp.Chem. 13(1992)371.
 */

GLOBAL_FUNCTION int
scf_charges_from_esp(centers,scf_info,charges,dipole,mulpts,
                     density,use_dip,outfile)
centers_t *centers;
scf_struct_t *scf_info;
double_vector_t *charges;
double_vector_t *dipole;
expts_t *mulpts;
double density;
int use_dip;
FILE *outfile;
{
  int errcod;
  int nat=centers->n;
  double_matrix_t A;
  double_vector_t B;


  errcod=allocbn_double_matrix(&A,"n1 n2",nat+4,nat+4);
  if(errcod!=0) {
    fprintf(outfile,"scf_charges_from_esp: could not alloc A matrix\n");
    return(-1);
    }
  zero_double_matrix(&A);

  errcod=allocbn_double_vector(&B,"n",nat+4);
  if(errcod!=0) {
    fprintf(outfile,"scf_charges_from_esp: could not alloc B vector\n");
    free_double_matrix(&A);
    return(-1);
    }
  zero_double_vector(&B);


 /* B is a vector of the form (b1,b2,...,bn,qtot,mu_x,mu_y,mu_z), so place
  * dipole components in it now.
  */

  B.d[nat+1]=dipole->d[0];
  B.d[nat+2]=dipole->d[1];
  B.d[nat+3]=dipole->d[2];

  errcod=atomic_charges_from_esp(charges,mulpts,&A,&B,density,use_dip,outfile);
  if(errcod!=0) {
    fprintf(outfile,"scf_charges_from_esp: atomic_charges_from_esp failed\n");
    free_double_vector(&B);
    free_double_matrix(&A);
    return(-1);
    }

  free_double_matrix(&A);
  free_double_vector(&B);
  return(0);
  }


LOCAL_FUNCTION int
atomic_charges_from_esp(charges,mulpts,A,B,density,use_dip,outfile)
double_vector_t *charges;
expts_t *mulpts;
double_matrix_t *A;
double_vector_t *B;
double density;
int use_dip;
FILE *outfile;
{
  int i,j,k,l,m,n;
  int errcod;
  int nat=B->n-4;
  int natri=nat*(nat+1)/2;
  int index=0;
  int me=mynode0();
  int nproc=numnodes0();
  int num_shell=4;           /* the number of surfaces to use */
  int npoint,mpoint;

  double pts_per_ang=density*0.529177060;
  double debye=2.54176548;
  double qtot=0.0;
  double Vi;
  double dist;
  double r_i,r_il,phi,theta;
  double rist,rict;
  double xp,yp,zp;
  double vdw_scale=1.4;
  double dx=0.0,dy=0.0,dz=0.0;
  double dtot;
  expts_t center;
  double_vector_t scr;

  center.n=nat;
  center.p=&mulpts->p[natri];

  errcod=allocbn_double_vector(&scr,"n",nat);
  if(errcod!=0) return(-1);

 /* determine total charge of the molecule and place in B */

  for(i=0; i<natri+nat; i++) qtot+=mulpts->p[i].charge;

  B->d[nat]=qtot;


 /* now, generate a large number of points 
  * use algorithm from Tasi, et al. for now
  *
  * we'll have 4 surfaces equal to 1.4, 1.6, 1.8, and 2.0 times the VDW radius
  */

  index=0;
  for(m=0; m < num_shell; m++) {
    for(i=0; i < nat; i++) {

  /* generate points on sphere around center i of radius scale*vdw_radius */
      r_i = vdw_scale*vdw_radius(center.p[i].charge);

      mpoint = (int) M_PI*r_i*pts_per_ang;
      if(mpoint%2) mpoint++;

      for(j=0; j <= mpoint; j++) {
        theta = M_PI*((double) j/(double)mpoint);

        rist=r_i*sin(theta);
        rict=r_i*cos(theta);

        npoint = (int) 2.0*M_PI*rist*pts_per_ang;
        if(!npoint || npoint%2) npoint++;

        for(k=0; k < npoint; k++,index++) {

          if(index%nproc!=me) continue;

          phi=2.0*M_PI*((double)k/(double)npoint);
          xp=rist*cos(phi)+center.p[i].r[0];
          yp=rist*sin(phi)+center.p[i].r[1];
          zp=rict+center.p[i].r[2];

       /* dist is the distance from point p to center l
        * if dist is within the sphere around center l, omit it
        */
          for(l=0; l < nat; l++) {
            if(l!=i) {
              dist=(xp-center.p[l].r[0])*(xp-center.p[l].r[0])+
                   (yp-center.p[l].r[1])*(yp-center.p[l].r[1])+
                   (zp-center.p[l].r[2])*(zp-center.p[l].r[2]);
              dist=sqrt(dist);
              r_il=vdw_scale*vdw_radius(center.p[l].charge);
              if(dist<r_il) break;
              }
            }
          if(l!=nat) continue;

        /* ok, this is a point on the surface.  Calculate the classical
         * ESP at this point, and then sum the contribution from this
         * point into A and B
         */

          Vi = esp(mulpts,xp,yp,zp);
          sum_A_B(A,B,Vi,&center,xp,yp,zp);
          }
        }
      }
    vdw_scale+=0.2;
    }

 /* sum contributions to A and B */

  gop1(B->d,nat,scr.d,'+',mtype_get());
  for(i=0; i < nat; i++) gop1(A->d[i],nat,scr.d,'+',mtype_get());

  free_double_vector(&scr);

 /* finish A matrix, so set of linear equations is:
  *
  * | A11 A12 ... A1n 1 x1 y1 z1 | | q1 |    | B1 |
  * |  .   .       .  .  . .  .  | | q2 |    | B2 |
  * |  .   .       .  .  . .  .  | | .  |    |  . |
  * |  .   .       .  .  . .  .  | | .  |    |  . |
  * | An1 An2 ... Ann 1 xn yn zn | | qn | =  | Bn |
  * |  1   1  ...  1  0  0  0  0 | |lamc|    |qtot|
  * |  x1  x2 ...  xn 0  0  0  0 | |lamx|    |mu_x|
  * |  y1  y2 ...  yn 0  0  0  0 | |lamy|    |mu_y|
  * |  z1  z2 ...  zn 0  0  0  0 | |lamz|    |mu_z|
  */

  for(i=0; i < nat; i++) {
    A->d[i][nat]=A->d[nat][i]=1.0;

    A->d[i][nat+1]=A->d[nat+1][i]=center.p[i].r[0];
    A->d[i][nat+2]=A->d[nat+2][i]=center.p[i].r[1];
    A->d[i][nat+3]=A->d[nat+3][i]=center.p[i].r[2];
    }

  if(use_dip) math_lin(A,B,B->n,1);
  else math_lin(A,B,B->n-3,1);

 /* copy charges into "charges" */

  for(i=0; i < charges->n; i++) charges->d[i]=B->d[i];

 /* re-calculate dipole moment from point charges as a test */

  for(i=0; i < nat; i++) {
    dx += B->d[i]*center.p[i].r[0];
    dy += B->d[i]*center.p[i].r[1];
    dz += B->d[i]*center.p[i].r[2];
    }

  dx*=debye;
  dy*=debye;
  dz*=debye;
  dtot = sqrt(dx*dx+dy*dy+dz*dz);

  if(me==0) {
    fprintf(outfile,"\n  dipole moment from ESP\n");
    fprintf(outfile,"     dx = %lf\n",dx);
    fprintf(outfile,"     dy = %lf\n",dy);
    fprintf(outfile,"     dz = %lf\n",dz);
    fprintf(outfile,"    tot = %lf\n",dtot);
    }

  return(0);
  }

LOCAL_FUNCTION double
vdw_radius(charge)
double charge;
{
  /* return van der Waals radius in a.u. */
  return 1.889726664*VDW_Radius[(int)charge];
  }

LOCAL_FUNCTION double
esp(mulpts,xp,yp,zp)
expts_t *mulpts;
double xp;
double yp;
double zp;
{
  int i;
  double rij;
  double *r;
  double v=0.0;

 /* the classical ESP = sum_i qi/|ri-rp| */

  for(i=0; i < mulpts->n; i++) {
    r=mulpts->p[i].r;
    rij = (xp-r[0])*(xp-r[0])+(yp-r[1])*(yp-r[1])+(zp-r[2])*(zp-r[2]);
    rij = 1.0/sqrt(rij);
    v+=rij*mulpts->p[i].charge;
    }
  return(v);
  }

LOCAL_FUNCTION VOID
sum_A_B(A,B,Vi,center,xp,yp,zp)
double_matrix_t *A;
double_vector_t *B;
double Vi;
expts_t *center;
double xp;
double yp;
double zp;
{
  int i,j,k;
  double rij,rik;
  double *r;

 /* Bk  = sum_i Vi/rik
  * Ajk = sum_i 1/(rij*rik)
  */

  for(k=0; k < center->n; k++) {
    r=center->p[k].r;
    rik=(xp-r[0])*(xp-r[0])+(yp-r[1])*(yp-r[1])+(zp-r[2])*(zp-r[2]);
    rik=1.0/sqrt(rik);

    B->d[k]+=Vi*rik;

    for(j=0; j <= k; j++) {
      r=center->p[j].r;
      rij=(xp-r[0])*(xp-r[0])+(yp-r[1])*(yp-r[1])+(zp-r[2])*(zp-r[2]);
      
      A->d[j][k]=A->d[k][j]+=rik/sqrt(rij); 
      }
    }
  }
