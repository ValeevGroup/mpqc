
#include <math.h>

#include <math/symmetry/pointgrp.h>
#include <util/misc/formio.h>
#include <iomanip.h>

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
SymmetryOperation::print(ostream& os) const
{
#if HAVE_IOS_FMTFLAGS
  ios::fmtflags oldf;
#else
  long oldf;
#endif
  oldf = os.setf(ios::left);

  os.setf(ios::fixed,ios::floatfield);
  os.setf(ios::right,ios::adjustfield);
  
  os << indent << "        1          2          3\n";
  os << indent << "  1  "
     << setw(10) << setprecision(7) << d[0][0] << " "
     << setw(10) << setprecision(7) << d[0][1] << " "
     << setw(10) << setprecision(7) << d[0][2] << "\n";
  os << indent << "  2  "
     << setw(10) << setprecision(7) << d[1][0] << " "
     << setw(10) << setprecision(7) << d[1][1] << " "
     << setw(10) << setprecision(7) << d[1][2] << "\n";
  os << indent << "  3  "
     << setw(10) << setprecision(7) << d[2][0] << " "
     << setw(10) << setprecision(7) << d[2][1] << " "
     << setw(10) << setprecision(7) << d[2][2] << "\n";

  os << endl;

  os.setf(oldf);
}

