
#include <util/class/class.h>
#include <util/state/state.h>

main(int argc, char *argv[])
{
  if(!argv[1]) exit(1);

  printf("avs output\n");
  printf("your comment here\n");
  printf("mopac stuff here\n");

  StateInBinXDR si(argv[1]);

  int natom;
  si.get(natom);

  printf("%5d%12.6f%12.6f%12.6f\n",natom, 0.0,      0.0,      0.0);
  printf("%5d%12.6f%12.6f%12.6f\n",natom, 0.0,      0.0,      0.755905);
  printf("%5d%12.6f%12.6f%12.6f\n",1,     0.0,      0.755905, 0.0);
  printf("%5d%12.6f%12.6f%12.6f\n",1,     0.755905, 0.0,      0.0);

  double x,y,z,Z;
  int n;

  for (int i=0; i < natom; i++) {
    si.get(n); si.get(Z); si.get(x); si.get(y); si.get(z);
    printf("%11d%17.7e%17.7e%17.7e%17.7e\n",n,Z,x,y,z);
  }

  si.get(Z); si.get(x); si.get(y); si.get(z);

  while(fabs(Z+99999.0) > 1.0e-3) {
    printf("%17.7e%17.7e%17.7e%17.7e\n",Z,x,y,z);
    si.get(Z); si.get(x); si.get(y); si.get(z);
  }

  exit(0);
}
