
#include <math.h>

#include <math/symmetry/pointgrp.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

////////////////////////////////////////////////////////////////////////

SymmetryOperation::SymmetryOperation()
{
  zero();
}

SymmetryOperation::~SymmetryOperation()
{
}

SymmetryOperation
SymmetryOperation::operate(const SymmetryOperation& r) const
{
  SymmetryOperation ret;
  for (int i=0; i < 3; i++)
    for (int j=0; j < 3; j++) {
      double t=0;
      for (int k=0; k < 3; k++)
        t += r.d[i][k]*d[k][j];
      ret.d[i][j] = t;
    }
  return ret;
}

SymmetryOperation
SymmetryOperation::sim_transform(const SymmetryOperation& r) const
{
  int i,j,k;
  SymmetryOperation ret,foo;
  
  // foo = r * d
  for (i=0; i < 3; i++) {
    for (j=0; j < 3; j++) {
      double t=0;
      for (k=0; k < 3; k++)
        t += r.d[i][k] * d[k][j];
      foo.d[i][j] = t;
    }
  }

  // ret = (r*d)*r~ = foo*r~
  for (i=0; i < 3; i++) {
    for (j=0; j < 3; j++) {
      double t=0;
      for (k=0; k < 3; k++)
        t += foo.d[i][k]*r.d[j][k];
      ret.d[i][j]=t;
    }
  }

  return ret;
}

// Clockwise rotation by 2pi/n degrees
void
SymmetryOperation::rotation(int n)
{
  double theta = (n) ? 2.0*M_PI/n : 2.0*M_PI;
  rotation(theta);
}

// Clockwise rotation by theta degrees
void
SymmetryOperation::rotation(double theta)
{
  zero();

  double ctheta = cos(theta);
  double stheta = sin(theta);

  d[0][0] = ctheta;
  d[0][1] = stheta;
  d[1][0] = -stheta;
  d[1][1] = ctheta;
  d[2][2] = 1.0;
}

void
SymmetryOperation::print(FILE* outfile) const
{
  fprintf(outfile,"        1          2          3\n");
  fprintf(outfile,"  1  %10.7f %10.7f %10.7f\n",d[0][0],d[0][1],d[0][2]);
  fprintf(outfile,"  2  %10.7f %10.7f %10.7f\n",d[1][0],d[1][1],d[1][2]);
  fprintf(outfile,"  3  %10.7f %10.7f %10.7f\n",d[2][0],d[2][1],d[2][2]);
  fprintf(outfile,"\n");
}

