
typedef int dmt_matrix;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <tmpl.h>

extern "C" {
#include <comm/picl/picl.h>
}

#include <util/class/class.h>
#include <util/state/state.h>
#include <util/keyval/keyval.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/simple.h>
#include <chemistry/molecule/symm.h>
#include <chemistry/molecule/simpleQCList.h>
#include <chemistry/molecule/symmQCList.h>
#include <math/nihmatrix/nihmatrix.h>
#include <math/nihmatrix/lmath.h>

extern "C" {
#include <chemistry/qc/dmtqc/libdmtqc.h>
#include <util/misc/libmisc.h>
#include <util/bio/libbio.h>
}

#include "mpqc_int.h"

extern "C" {
 int mynode0();
 int numnodes0();
 void gop1(double*,int,double*,char,int);
}

static double esp(expts_t*, double , double , double);
static double espr(centers_t *, expts_t*, dmt_matrix, double , double , double );
static void sum_A_B(DMatrix&, DVector&, double, expts_t*,double,double,double);
static int atomic_charges_from_esp(centers_t *centers, Molecule& mol,
  dmt_matrix,
  double_vector_t *charges, expts_t *mulpts,
  DMatrix& A, DVector& B,
  double density, int use_dip, FILE *outfile, RefKeyVal keyval);

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

int
Scf_charges_from_esp(centers_t *centers, dmt_matrix scf_vec,
    double_vector_t *charges, double_vector_t *dipole, expts_t *mulpts,
    double density, int use_dip, FILE *outfile, RefKeyVal keyval)
{

  int errcod;
  int nat=centers->n;

#define USE_DIPOLE 0
#if USE_DIPOLE
  DMatrix A(nat+4,nat+4);
  DVector B(nat+4);
#else
  DMatrix A(nat+1,nat+1);
  DVector B(nat+1);
#endif

  A.zero(); B.zero();

 /* B is a vector of the form (b1,b2,...,bn,qtot,mu_x,mu_y,mu_z), so place
  * dipole components in it now.
  */

#if USE_DIPOLE
  B[nat+1]=dipole->d[0];
  B[nat+2]=dipole->d[1];
  B[nat+3]=dipole->d[2];
#endif

 // now lets make another molecule

  Molecule mol;

  double charge=0;
  int natr=nat*(nat+1)/2;
  for(int i=0; i < nat; i++) {
    centers->center[i].r[0] = mulpts->p[natr+i].r[0];
    centers->center[i].r[1] = mulpts->p[natr+i].r[1];
    centers->center[i].r[2] = mulpts->p[natr+i].r[2];

    AtomicCenter ac(centers->center[i].atom,
		    centers->center[i].r[0],
		    centers->center[i].r[1],
		    centers->center[i].r[2]);
    mol.add_atom(i,ac);
    charge += centers->center[i].charge;
    }

  int ndoc = (int) charge / 2;

  int_initialize_1e(0,0,centers,centers);
  int_initialize_offsets1(centers,centers);

  dmt_matrix pmat = dmt_create("density",centers->nfunc,SCATTERED);
  dmt_density(scf_vec,ndoc,pmat);
  dmt_scale(pmat,2.0);

  atomic_charges_from_esp(centers,mol,pmat,charges,mulpts,A,B,
                          density,use_dip,outfile,keyval);

  int_done_offsets1(centers,centers);
  int_done_1e();

  dmt_free(pmat);

  return(0);
  }


static int
atomic_charges_from_esp(centers_t *centers, Molecule& mol, dmt_matrix pmat,
  double_vector_t *charges, expts_t *mulpts,
  DMatrix& A, DVector& B,
  double density, int use_dip, FILE *outfile, RefKeyVal keyval)
{
  int i,j,k,l,m,n;
  int errcod;
  int nat=centers->n;
  int natri=nat*(nat+1)/2;
  int index=0;
  int me=mynode0();
  int nproc=numnodes0();
  int num_shell=4;           /* the number of surfaces to use */
  int npoint,mpoint;

  double bohr = 0.52917706;
  double pts_per_ang=density*0.529177060;
  double debye=2.54176548;
  double qtot=0.0;
  double Vi;
  double dist;
  double r_i,r_il,phi,theta;
  double rist,rict;
  double xp,yp,zp;
  double vdw_scale=1.4;
  double vdw_spacing=0.2;
  double dx=0.0,dy=0.0,dz=0.0;
  double dtot;
  expts_t center;
  double_vector_t scr;

  center.n=nat;
  center.p=&mulpts->p[natri];

  errcod=allocbn_double_vector(&scr,"n",nat);
  if(errcod!=0) return(-1);

 /* read in some input junk */
  if (me==0) {
    if (keyval->exists("num_shell"))
      num_shell=keyval->intvalue("num_shell");

    if (keyval->exists("vdw_scale"))
      vdw_scale = keyval->doublevalue("vdw_scale");

    if (keyval->exists("vdw_spacing"))
      vdw_spacing = keyval->doublevalue("vdw_spacing");
    }

  bcast0(&num_shell,sizeof(int),mtype_get(),0);
  bcast0(&vdw_scale,sizeof(double),mtype_get(),0);
  bcast0(&vdw_spacing,sizeof(double),mtype_get(),0);

 /* set up file for avs */

  //StateOutBin avsfil((FILE*)0);
  StateOutBinXDR avsfil((FILE*)0);

  if (me==0) {
    avsfil.open("esp.dat","w");

    avsfil.put(mol.natom());

    for (i=0; i < mol.natom(); i++) {
      avsfil.put(mol[i].element().number());
      avsfil.put(mol[i].element().charge());
      avsfil.put(center.p[i].r[0]);
      avsfil.put(center.p[i].r[1]);
      avsfil.put(center.p[i].r[2]);
      }
    avsfil.flush();
    }

 /* determine total charge of the molecule and place in B */

  for(i=0; i<natri+nat; i++) qtot+=mulpts->p[i].charge;

  B[nat]=qtot;

 /* now, generate a large number of points 
  * use algorithm from Tasi, et al. for now
  *
  * we'll have 4 surfaces equal to 1.4, 1.6, 1.8, and 2.0 times the VDW radius
  */

  index=0; int npoints=0;
  for(m=0; m < num_shell; m++) {
    for(i=0; i < nat; i++) {

  /* generate points on sphere around center i of radius scale*vdw_radius */
      r_i = vdw_scale*mol[i].element().vdw_radius();

      mpoint = (int) (M_PI*r_i*pts_per_ang);
      if(mpoint%2) mpoint++;

      for(j=0; j <= mpoint; j++) {
        theta = M_PI*((double) j/(double)mpoint);

        rist=r_i*sin(theta);
        rict=r_i*cos(theta);

        npoint = (int) (2.0*M_PI*rist*pts_per_ang);
        //if(!npoint || npoint%2) npoint++; // makes a pretty picture
        if (!npoint) npoint++;

        for(k=0; k < npoint; k++,index++) {

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
              r_il=vdw_scale*mol[l].element().vdw_radius();
              if(dist<r_il) break;
              }
            }
          if(l!=nat) continue;

        /* ok, this is a point on the surface.  Calculate the classical
         * ESP at this point, and then sum the contribution from this
         * point into A and B
         */

#define CLASSICAL 0
#if CLASSICAL
	  Vi = esp(mulpts,xp,yp,zp);

          if(index%nproc==me) {
	    sum_A_B(A,B,Vi,&center,xp,yp,zp);
            }

          if (me==0) {
	    avsfil.put(Vi);
            avsfil.put(xp);
            avsfil.put(yp);
            avsfil.put(zp);
            avsfil.flush();
            }
#else
	  Vi = espr(centers,&center,pmat,xp,yp,zp);
	  sum_A_B(A,B,Vi,&center,xp,yp,zp);
          if (me==0) {
	    avsfil.put(Vi);
            avsfil.put(xp);
            avsfil.put(yp);
            avsfil.put(zp);
            avsfil.flush();
            }
          sync0();
#endif
          npoints++;
          }
        }
      }
    vdw_scale+=vdw_spacing;
    }

  if (me==0) {
    vdw_scale=-99999.0;
    avsfil.put(vdw_scale);
    avsfil.put(vdw_scale);
    avsfil.put(vdw_scale);
    avsfil.put(vdw_scale);
    }

 /* sum contributions to A and B */

#if CLASSICAL
  gop1(B.pointer(),nat,scr.d,'+',mtype_get());
  for (i=0; i < nat; i++) gop1(A[i],nat,scr.d,'+',mtype_get());
#endif

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
    A(i,nat)=A(nat,i)=1.0;

#if USE_DIPOLE
    A(i,nat+1) = A(nat+1,i) = center.p[i].r[0];
    A(i,nat+2) = A(nat+2,i) = center.p[i].r[1];
    A(i,nat+3) = A(nat+3,i) = center.p[i].r[2];
#endif
    }
 
  A.solve_lin(B);

 /* copy charges into "charges" */

  for(i=0; i < charges->n; i++) charges->d[i]=B[i];

 /* re-calculate dipole moment from point charges as a test */

  for(i=0; i < nat; i++) {
    dx += B[i]*center.p[i].r[0];
    dy += B[i]*center.p[i].r[1];
    dz += B[i]*center.p[i].r[2];
    }

  dx*=debye;
  dy*=debye;
  dz*=debye;
  dtot = sqrt(dx*dx+dy*dy+dz*dz);

  if(me==0) {
    fprintf(outfile,"\n  %d ESP points generated\n",npoints);
    fprintf(outfile,"\n  dipole moment from ESP\n");
    fprintf(outfile,"     dx = %lf\n",dx);
    fprintf(outfile,"     dy = %lf\n",dy);
    fprintf(outfile,"     dz = %lf\n",dz);
    fprintf(outfile,"    tot = %lf\n",dtot);
    }

  return(0);
  }



static double
esp(expts_t *mulpts, double xp, double yp, double zp)
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

static void
sum_A_B(DMatrix& A, DVector& B, double Vi, expts_t *center,
        double xp, double yp, double zp)
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

    B[k]+=Vi*rik;

    for(j=0; j <= k; j++) {
      r=center->p[j].r;
      rij=(xp-r[0])*(xp-r[0])+(yp-r[1])*(yp-r[1])+(zp-r[2])*(zp-r[2]);
      
      A[j][k]=A[k][j]+=rik/sqrt(rij); 
      }
    }
  }


static double
espr(centers_t *centers, expts_t *center,
    dmt_matrix pmat, double xp, double yp, double zp)
{
  int i;
  double rij;
  double *r;
  double vc=0,vq=0.0;

  double pos[3];
  double buff[36];
  pos[0]=xp; pos[1]=yp; pos[2]=zp;

 /* the classical ESP = sum_i qi/|ri-rp| */

  for(i=0; i < centers->n; i++) {
   // r=centers->center[i].r;
    r=center->p[i].r;
    rij = (xp-r[0])*(xp-r[0])+(yp-r[1])*(yp-r[1])+(zp-r[2])*(zp-r[2]);
    rij = 1.0/sqrt(rij);
    vc += rij*centers->center[i].charge;
    }

  for(int nl=0; nl < dmt_nlocal(pmat); nl++) {
    int ib,jb,isz,jsz,ist,jst;
    double *lblk;

    dmt_get_block(pmat,nl,&ib,&jb,&lblk);
    dmt_describe_block(pmat,ib,&ist,&isz);
    dmt_describe_block(pmat,jb,&jst,&jsz);

    double *ppos=pos;
    double one=1;
    int_shell_point_charge(centers,centers,buff,ib,jb,1,&one,&ppos);

    double scal = (ib==jb) ? 1.0 : 2.0 ;

    for(i=0; i < isz*jsz; i++) vq += scal*lblk[i]*buff[i];
    }

  gsum0(&vq,1,5,mtype_get(),0);

  return vc+vq;
  }
