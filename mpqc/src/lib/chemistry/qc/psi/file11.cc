extern "C" {
#include<stdio.h>
int read_file11(int natom, int readto,
      char *label, char *theory, char *dertype,
      double *energy, double *X, double *Y, double *Z, int *AN,
      double **grad, int *ngrad);
void print_file11(char *label, int natom, double energy, double *X,
      double *Y, double *Z, int *AN, FILE *fpo);
}
#include "file11.h"

FILE11::FILE11(int num)
{
  char junk[132];

  FILE *fp = fopen("file11.dat","r");
  fgets(junk, 132, fp);
  fscanf(fp, "%d", &nat);
  fclose(fp);

  for(int i = 0; i<3; i++){
     coord_[i] = new double[nat];
     grad_[i] = new double[nat];
     }
 charges_ = new int[nat];
  
 int errcod = read(num);
 if (!errcod)
   fprintf(stderr, "failed to read gradient #%d from file11\n", num);
}
FILE11::~FILE11()
{
  for(int i = 0; i<3; i++){
    delete[] coord_[i];
    delete[] grad_[i];
    }
  delete[] charges_;
}
int FILE11::read(int num)
{
  return read_file11(nat, num, label_, theory_, dertype_, &energy_, 
                coord_[0], coord_[1], coord_[2], charges_, grad_, &ngrad_);
}
void FILE11::print()
{
  print_file11(label_, nat, energy_, coord_[0], coord_[1], coord_[2], 
      charges_, stdout);
}
