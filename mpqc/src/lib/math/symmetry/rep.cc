
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <math/symmetry/pointgrp.h>

/////////////////////////////////////////////////////////////////////////

SymRep::SymRep(int i) :
  n(i)
{
  zero();
}

SymRep::SymRep(const SymmetryOperation& so) :
  n(3)
{
  memset(d,0,sizeof(double)*25);
  for (int i=0; i < 3; i++)
    for (int j=0; j < 3; j++)
      d[i][j] = so[i][j];
}

SymRep::~SymRep()
{
  n=0;
}

SymRep::operator SymmetryOperation() const
{
  if (n != 3) {
    fprintf(stderr,"SymRep::operator SymmetryOperation(): ");
    fprintf(stderr,"trying to cast to symop when n == %d\n",n);
    abort();
  }
    
  SymmetryOperation so;

  for (int i=0; i < 3; i++)
    for (int j=0; j < 3; j++)
      so[i][j] = d[i][j];

  return so;
}

SymRep
SymRep::operate(const SymRep& r) const
{
  if (r.n != n) {
    fprintf(stderr,"SymRep::operate(): dimensions don't match\n");
    fprintf(stderr,"  %d != %d\n",r.n,n);
    abort();
  }
  
  SymRep ret(n);

  for (int i=0; i < n; i++) {
    for (int j=0; j < n; j++) {
      double t=0;
      for (int k=0; k < n; k++)
        t += r[i][k] * d[k][j];
      ret[i][j] = t;
    }
  }

  return ret;
}

SymRep
SymRep::sim_transform(const SymRep& r) const
{
  int i,j,k;

  if (r.n != n) {
    fprintf(stderr,"SymRep::sim_transform(): dimensions don't match\n");
    fprintf(stderr,"  %d != %d\n",r.n,n);
    abort();
  }
  
  SymRep ret(n), foo(n);

  // foo = r * d
  for (i=0; i < n; i++) {
    for (j=0; j < n; j++) {
      double t=0;
      for (k=0; k < n; k++)
        t += r[i][k] * d[k][j];
      foo[i][j] = t;
    }
  }

  // ret = (r*d)*r~ = foo*r~
  for (i=0; i < n; i++) {
    for (j=0; j < n; j++) {
      double t=0;
      for (k=0; k < n; k++)
        t += foo[i][k]*r[j][k];
      ret[i][j]=t;
    }
  }

  return ret;
}
  
void
SymRep::rotation(int nt)
{
  double theta = (nt) ? 2.0*M_PI/nt : 2.0*M_PI;
  rotation(theta);
}

void
SymRep::rotation(double theta)
{
  zero();

  double ctheta = cos(theta);
  double stheta = sin(theta);
  double c2theta = cos(2*theta);
  double s2theta = sin(2*theta);

  switch (n) {
  case 1:
    d[0][0] = 1.0;
    break;
    
  case 3:
    d[2][2] = 1.0;
    // deliberately fall through

  case 4:
  case 2:
    d[0][0] = ctheta;
    d[0][1] = stheta;
    d[1][0] = -stheta;
    d[1][1] = ctheta;
    
    // this is ok since d is hardwired
    d[2][2] = c2theta;
    d[2][3] = -s2theta;
    d[3][2] = s2theta;
    d[3][3] = c2theta;
    break;
    
  case 5:
    d[0][0] = 1.0;
    
    d[1][1] = c2theta;
    d[1][2] = s2theta;
    d[2][1] = -s2theta;
    d[2][2] = c2theta;

    d[3][3] = ctheta;
    d[3][4] = -stheta;
    d[4][3] = stheta;
    d[4][4] = ctheta;
    break;

  default:
    fprintf(stderr,"SymRep::rotation(): n > 5 (%d)\n",n);
    abort();
  }
  
}

void
SymRep::print(FILE* outfile) const
{
  int i;
  
  for (i=0; i < n; i++)
    fprintf(outfile,"%11d",i+1);
  fprintf(outfile,"\n");
  
  for (i=0; i < n; i++) {
    fprintf(outfile,"%3d ",i+1);
    for (int j=0; j < n; j++)
      fprintf(outfile," %10.7f",d[i][j]);
    fprintf(outfile,"\n");
  }
  fprintf(outfile,"\n");
}
