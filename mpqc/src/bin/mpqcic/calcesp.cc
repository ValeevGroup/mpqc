
#include <util/class/class.h>
#include <util/state/state.h>
#include <math/nihmatrix/nihmatrix.h>
#include <chemistry/qc/dmtscf/masses.h>

struct points {
  double x;
  double y;
  double z;
  double Z;
} pts,*atoms,*mulpts;

DMatrix transmat(3,3);
DVector com(3);

static void readinp(FILE* inp, int& natom, int& npts);
static int read_mul(char *path);
static void check_geom(int,int);
static void center_of_mass(struct points *at, int natom);

int
main(int argc, char *argv[])
{
  if(!argv[1]) exit(1);

  FILE *inp = fopen(argv[1],"r");
  if (!inp) exit(2);

 // if second argument given, read in point charges
  int nmul=0;
  if (argc==3) nmul = read_mul(argv[2]);

 // read in some header stuff from inp
  int natom=0,npts=0;
  readinp(inp,natom,npts);

  DMatrix A(natom+1,natom+1);
  DVector B(natom+1);

 // read in geometry
  int an;
  for (int i=0; i < natom; i++)
    fscanf(inp,"%d%lf%lf%lf%lf\n",&an,&atoms[i].Z,
                                      &atoms[i].x,
                                      &atoms[i].y,
                                      &atoms[i].z);

 // make sure geometry in atoms is same as in mulpts
  check_geom(natom,nmul);


 // and then make A matrix and B vector
  A.zero(); B.zero();

  for (i=0; i < npts; i++) {
    fscanf(inp,"%lf%lf%lf%lf\n",&pts.Z,&pts.x,&pts.y,&pts.z);

    if (mulpts) {
      pts.Z=0;
      double x = pts.x - com[0];
      double y = pts.y - com[1];
      double z = pts.z - com[2];

      pts.x = transmat(0,0)*x + transmat(1,0)*y + transmat(2,0)*z;
      pts.y = transmat(0,1)*x + transmat(1,1)*y + transmat(2,1)*z;
      pts.z = transmat(0,2)*x + transmat(1,2)*y + transmat(2,2)*z;

      for (int j=0; j < nmul; j++) {
	double rij = sqrt((mulpts[j].x-pts.x)*(mulpts[j].x-pts.x) +
			  (mulpts[j].y-pts.y)*(mulpts[j].y-pts.y) +
			  (mulpts[j].z-pts.z)*(mulpts[j].z-pts.z));
        pts.Z += mulpts[j].Z/rij;
        }
      }

    for (int k=0; k < natom; k++) {
      double rik = sqrt((atoms[k].x-pts.x)*(atoms[k].x-pts.x) +
                        (atoms[k].y-pts.y)*(atoms[k].y-pts.y) +
                        (atoms[k].z-pts.z)*(atoms[k].z-pts.z));

      B[k] += pts.Z/rik;

      for (int j=0; j <= k; j++) {
	double rij = sqrt((atoms[j].x-pts.x)*(atoms[j].x-pts.x) +
			  (atoms[j].y-pts.y)*(atoms[j].y-pts.y) +
			  (atoms[j].z-pts.z)*(atoms[j].z-pts.z));

        A[k][j] = A[j][k] += 1.0/(rij*rik);
        }
      }
    }

  for(i=0; i < natom; i++) A(i,natom)=A(natom,i)=1.0;

  A.solve_lin(B);

  printf (" charges\n");

  for (i=0; i < natom; i++) 
    printf(" %10d %20.10f\n",i+1,B[i]);

  exit(0);
}


static int
read_mul(char *path)
{
  FILE *inp2 = fopen(path,"r");
  if (!inp2) exit(3);

  int nmul;
  fscanf(inp2,"%d\n",&nmul);

  mulpts = new struct points[nmul];

  for (int i=0; i < nmul; i++) 
    fscanf(inp2,"%lf%lf%lf%lf\n",&mulpts[i].Z,&mulpts[i].x,&mulpts[i].y,
                                   &mulpts[i].z);
  fclose(inp2);

  return nmul;
}

static void
readinp(FILE* inp, int& natom, int& npts)
{
  char foo[80];
  
  fgets(foo,79,inp);
  fgets(foo,79,inp);
  fgets(foo,79,inp);

  int an;
  double x,y,z;

  fscanf(inp,"%d%lf%lf%lf\n",&natom,&x,&y,&z);

  fscanf(inp,"%d%lf%lf%lf\n",&npts,&x,&y,&z);
  fscanf(inp,"%d%lf%lf%lf\n",&an,&x,&y,&z);
  fscanf(inp,"%d%lf%lf%lf\n",&an,&x,&y,&z);

  atoms = new struct points[natom];
}

static void
check_geom(int natom, int nmul)
{
  int i;
  int swapx=0,swapy=0,swapz=0;

  // if we're not recalculating esp, then don't need to do any of this
  if (!mulpts) return;

  // first let's make sure mulpts geometry is c.o.m.

  struct points *at = &mulpts[nmul-natom];

  center_of_mass(at,natom);

  for (i=0; i < nmul-natom; i++) {
    double x = mulpts[i].x - com[0];
    double y = mulpts[i].y - com[1];
    double z = mulpts[i].z - com[2];

    mulpts[i].x = transmat(0,0)*x + transmat(1,0)*y + transmat(2,0)*z;
    mulpts[i].y = transmat(0,1)*x + transmat(1,1)*y + transmat(2,1)*z;
    mulpts[i].z = transmat(0,2)*x + transmat(1,2)*y + transmat(2,2)*z;
  }

  center_of_mass(atoms,natom);

  for (i=0; i < natom; i ++) {
    if ((fabs(at[i].x) > 1.0e-8) && fabs(at[i].x+atoms[i].x) < 0.1) swapx=1;
    if ((fabs(at[i].y) > 1.0e-8) && fabs(at[i].y+atoms[i].y) < 0.1) swapy=1;
    if ((fabs(at[i].z) > 1.0e-8) && fabs(at[i].z+atoms[i].z) < 0.1) swapz=1;
    }

  for (i=0; i < nmul; i ++) {
    if (swapx) mulpts[i].x *= -1;
    if (swapy) mulpts[i].y *= -1;
    if (swapz) mulpts[i].z *= -1;
    }

  for (i=0; i < natom; i++) {
    if (fabs(at[i].x-atoms[i].x) > 1.0e-1 ||
        fabs(at[i].y-atoms[i].y) > 1.0e-1 ||
        fabs(at[i].z-atoms[i].z) > 1.0e-1 ) {

      printf("poor fit %d\n",i+1);
      }
    }
}

static void
center_of_mass(struct points *at, int natom)
{
  int i;
  double X,Y,Z,M;
  double ausq_to_angsq = 0.2800283608302436;

  DVector x(natom);
  DVector y(natom);
  DVector z(natom);
  DVector m(natom);


  X=Y=Z=M=0.0;
  for(i=0; i < natom; i++) {
    x[i]=at[i].x;
    y[i]=at[i].y;
    z[i]=at[i].z;
    m[i]=Atomic_Mass[(int)at[i].Z];
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

  DVector evals(3);
  DMatrix inert(3,3);
  inert.zero();

  for (i=0; i < natom ; i++) {
    x[i] = (x[i]-X);
    y[i] = (y[i]-Y);
    z[i] = (z[i]-Z);

    inert(0,0) += ausq_to_angsq*m[i]*(y[i]*y[i] + z[i]*z[i]);
    inert(1,0) -= ausq_to_angsq*m[i]*x[i]*y[i];
    inert(1,1) += ausq_to_angsq*m[i]*(x[i]*x[i] + z[i]*z[i]);
    inert(2,0) -= ausq_to_angsq*m[i]*x[i]*z[i];
    inert(2,1) -= ausq_to_angsq*m[i]*y[i]*z[i];
    inert(2,2) += ausq_to_angsq*m[i]*(x[i]*x[i] + y[i]*y[i]);
    }

  inert(0,1) = inert(1,0);
  inert(0,2) = inert(2,0);
  inert(1,2) = inert(2,1);

  inert.print("I matrix");

  inert.diagonalize(evals,transmat);

  transmat.print("A matrix");

  evals.print("principal moments");

  for (i=0; i < natom ; i++) {
    at[i].x = transmat(0,0)*x[i] + transmat(1,0)*y[i] + transmat(2,0)*z[i];
    at[i].y = transmat(0,1)*x[i] + transmat(1,1)*y[i] + transmat(2,1)*z[i];
    at[i].z = transmat(0,2)*x[i] + transmat(1,2)*y[i] + transmat(2,2)*z[i];
  }
}
