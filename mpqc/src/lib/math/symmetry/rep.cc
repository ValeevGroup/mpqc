
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <math/symmetry/pointgrp.h>
#include <util/misc/formio.h>
#include <iomanip.h>

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
    cerr << indent << "SymRep::operator SymmetryOperation(): "
         << "trying to cast to symop when n == " << n << endl;
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
    cerr << indent << "SymRep::operate(): dimensions don't match: "
         << r.n << " != " << n << endl;
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
    cerr << indent << "SymRep::sim_transform(): dimensions don't match: "
         << r.n << " != " << n << endl;
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
SymRep::sigma_h()
{
  unit();

  if (n==3) {
    d[2][2] = -1.0;
  } else if (n==5) {
    d[3][3] = d[4][4] = -1.0;
  }
}
  
void
SymRep::sigma_xz()
{
  unit();

  if (n==2 || n==3 || n==4) {
    d[1][1] = -1.0;
    if (n==4)
      d[2][2] = -1.0;
  } else if (n==5) {
    d[2][2] = d[4][4] = -1.0;
  }
}

void
SymRep::sigma_yz()
{
  unit();

  if (n==2 || n==3 || n==4) {
    d[0][0] = -1.0;
    if (n==4)
      d[3][3] = -1.0;
  } else if (n==5) {
    d[2][2] = d[3][3] = -1.0;
  }
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
    d[0][0] = ctheta;
    d[0][1] = stheta;
    d[1][0] = -stheta;
    d[1][1] = ctheta;
    d[2][2] = 1.0;
    break;

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
    cerr << indent << "SymRep::rotation(): n > 5 (" << n << ")\n";
    abort();
  }
  
}

void
SymRep::c2_x()
{
  i();

  if (n==2 || n==3 || n==4) {
    d[0][0] = 1.0;
    if (n==4)
      d[3][3] = 1.0;
  } else if (n==5) {
    d[0][0] = d[1][1] = d[4][4] = 1.0;
  }
}
  
void
SymRep::c2_y()
{
  i();

  if (n==2 || n==3 || n==4) {
    d[1][1] = 1.0;
    if (n==4)
      d[2][2] = 1.0;
  } else if (n==5) {
    d[0][0] = d[1][1] = d[3][3] = 1.0;
  }
}
  

void
SymRep::print(ostream& os) const
{
  int i;

#ifdef HAVE_IOS_FMTFLAGS
  ios::fmtflags oldf;
#else
  long oldf;
#endif
  oldf = os.setf(ios::left);

  os.setf(ios::fixed,ios::floatfield);
  os.setf(ios::right,ios::adjustfield);
  
  os << indent;
  for (i=0; i < n; i++) os << setw(11) << i+1;
  os << endl;
  
  for (i=0; i < n; i++) {
    os << indent << setw(3) << i+1 << " ";
    for (int j=0; j < n; j++)
      os << " " << setw(10) << setprecision(10) << d[i][j];
    os << endl;
  }
  os << endl;

  os.setf(oldf);
}
