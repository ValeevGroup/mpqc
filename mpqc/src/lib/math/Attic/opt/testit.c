#include <math.h>
#define MAX_ATOMS 1000

/* this is the number of energy evaluations */
int n_energy_eval=0;
int n_energy_eval_tot=0;

/* define a 3d data type */
typedef struct {
  double x;
  double y;
  double z;
} coord_struct;

double calc_energy(coord_struct *, int);
void   angle_to_cart(double *, coord_struct *, int);

main()
{
  int n_atoms,i,j,n_trials,idum= -1;
  double angle[MAX_ATOMS];
  double ran3();
  double fitness, energy;
  coord_struct coords[MAX_ATOMS];
  scanf("%d %d",&n_atoms,&n_trials);
  
  for (j=0;j<n_trials;j++)
    {
      for (i=0;i<abs(n_atoms)-2;i++)
	{
	  angle[i]=360.*ran3(&idum);
	  printf("%4.0lf ",angle[i]);
	}
      printf("\n");
    }
  /* convert the angles to radians to avoid confusion later */
  for (i=0;i<n_atoms-2;i++) {
    angle[i]=angle[i]*M_PI/180.0;
    printf("%12.4e",angle[i]);
  }
  printf("\n");

  angle_to_cart(angle,coords,n_atoms);

  energy = calc_energy(coords,n_atoms);
  
  printf("energy=%lf\n",energy);
}

ljeval(int mode, int n, ColumnVector x, double& f, ColumnVector& g)
{
/* Compute f, grad */

}
