
#include <stdlib.h>
#include "nihmatrix.h"
#include <math/topology/point.h>
#include <util/keyval/keyval.h>

ARRAY_def(DVector);
DescribedClass_REF_def(DVector);

#define CLASSNAME DVector
#define PARENTS virtual public SavableState
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
DVector::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}

DVector::DVector() : n(0), d(0) {}

DVector::DVector(int sz) : n(sz)
{
  d=new double[n];
}

DVector::DVector(int sz, double *dp) : n(sz), d(dp) {}

DVector::DVector(const DVector& da) : n(0), d(0)
{
  *this=da;
}

DVector::DVector(RefPoint& da) : n(0), d(0)
{
  *this=da;
}

DVector::~DVector()
{
  n=0;
  if(d) delete[] d;
  d=0;
}

void DVector::resize(int nn)
{
  if(!nn) { init(); return; }

  if(n!=nn) {
    n=nn; if(d) delete[] d;
    d = new double[n];
    }
  if(!d)
    d = new double[n];
  }

DVector& DVector::operator=(const DVector& dv)
{
  resize(dv.n);
  for(int i=0; i < n; i++) d[i] = dv.d[i];
  return *this;
  }

DVector& DVector::operator=(RefPoint& dv)
{
  resize(dv->dimension());
  for(int i=0; i < n; i++) d[i] = dv->operator[](i);
  return *this;
  }

////////////////////////////////////////////////////////////////

DVector& DVector::operator+=(const double s)
{
  for(int i=0; i < n; i++) d[i] += s;
  return *this;
  }

DVector& DVector::operator*=(const double s)
{
  for(int i=0; i < n; i++) d[i] *= s;
  return *this;
  }

DVector& DVector::operator+=(const DVector& dv)
{
  if(dv.n!=n)
    err_quit("DVector::operator+=: vectors not same length %d %d",n,dv.n);

  for(int i=0; i < n; i++) d[i] += dv.d[i];
  return *this;
  }

DVector& DVector::operator-=(const DVector& dv)
{
  if(dv.n!=n)
    err_quit("DVector::operator-=: vectors not same length %d %d",n,dv.n);

  for(int i=0; i < n; i++) d[i] -= dv.d[i];
  return *this;
  }

DVector& DVector::operator*=(const DVector& dv)
{
  if(dv.n!=n)
    err_quit("DVector::operator*=: vectors not same length %d %d",n,dv.n);

  for(int i=0; i < n; i++) d[i] *= dv.d[i];
  return *this;
  }

/////////////////////////////////////////////////////////////////////

void
DVector::save_data_state(StateOut& so)
{
  so.put(n); so.put(d,n);
}

DVector::DVector(KeyVal&k)
{
  n = k.count();
  d = new double[n];
  for (int i=0; i<n; i++) {
      d[i] = k.doublevalue(i);
    }
}

DVector::DVector(StateIn&si):
  SavableState(si,class_desc_)
{
  si.get(n); si.get(d);
}

/////////////////////////////////////////////////////////////////////

const double DVector::maxval() const
{
  double t=0,tt;
  for(int i=0; i < n; i++) if((tt=fabs(d[i]))>t) t=tt;
  return t;
  }

double DVector::dot(const DVector& dv) const
{
  if(n!=dv.n)
    err_quit("DVector::dot(): dimensions not equal");

  double t=0; double *dp=d; double *dvp=dv.d;

  for(int i=0; i<n; i++) t += *dp++ * *dvp++;
  return t;
  }

double DVector::norm() const
{

  double t=0; double *dp=d;

  for(int i=0; i<n; i++,dp++) t += *dp * *dp;
  return sqrt(t);
  }

void DVector::normalize()
{

  double t=1.0/norm(); double *dp=d;

  for(int i=0; i<n; i++,dp++) *dp = *dp * t;
  }

#ifdef __GNUC__
DMatrix DVector::ccross(const DVector& dv) const return t(n,dv.n);
{
  for(int i=0; i < n; i++)
    for(int j=0; j < dv.n; j++)
      t(i,j) = d[i]*dv.d[j];
}
#else
DMatrix DVector::ccross(const DVector& dv) const
{
  DMatrix t(n,dv.n);

  for(int i=0; i < n; i++)
    for(int j=0; j < dv.n; j++)
      t(i,j) = d[i]*dv.d[j];

  return t;
  }
#endif

DVector DVector::cross(const DVector&v) const
#ifdef __GNUC__
     return result(v.n);
#endif
{
  if (v.n != 3 || n != 3) {
      fprintf(stderr,"DVector::cross: only handles dim 3\n");
      abort();
    }
#ifndef __GNUC__
  DVector result(v.n);
#endif
  result[0] = d[1]*v.d[2]-d[2]*v.d[1];
  result[1] = d[2]*v.d[0]-d[0]*v.d[2];
  result[2] = d[0]*v.d[1]-d[1]*v.d[0];
#ifndef __GNUC__
  return result;
#endif
}

void DVector::print(ostream& os) const { print(0,os); }

void DVector::print(const char *title, ostream& os, int prec) const
{
os << "in print n==" << this->n << ", d==" << d << endl;

  int ii,jj,kk,nn,ll;
  int i,j;
  int lwidth,width;
  double max=this->maxval();

os << "afger maxval\n";

  max=(max==0.0)?1.0:log10(max);
  if(max < 0.0) max=1.0;

  lwidth = prec+5+(int) max; width = 75/lwidth;

  os.setf(ios::fixed,ios::floatfield); 
  os.precision(prec);
  os.setf(ios::right,ios::adjustfield);

  if(title) os << "\n" << title << "\n";
  else os << "\n";

  if(n==0) { os << " empty vector\n"; return; }

  for(ii=jj=0;;) {
    ii++; jj++; kk=width*jj;
    nn=(n>kk)?kk:n;
    ll=2*(nn-ii+1)+1;

 // print column indices
    for(i=ii; i <= nn; i++) { os.width(lwidth); os << i; }
    os << "\n";

    for(j=ii-1; j < nn; j++) { os.width(lwidth); os << d[j]; }
    os << "\n\n";

    if(n<=kk) { os.flush(); return; }
    ii=kk;
    }
  }

////////////////////////////////////////////////////////////////////////

double dot(const DVector& v1, const DVector& v2)
{
  return v1.dot(v2);
  }

#ifdef __GNUC__
DVector operator+(const DVector& dv, const double s) return t;
{
  t = dv;
  t += s;
}

DVector operator+(double s, const DVector& dv) return t;
{
  t = dv;
  t += s;
}

DVector operator*(const DVector& dv, const double s) return t;
{
  t = dv;
  t *= s;
}

DVector operator*(double s, const DVector& dv) return t;
{
  t = dv;
  t *= s;
}

DVector operator+(const DVector& v1, const DVector& v2) return t;
{
  t = v1;
  t += v2;
}

DVector operator-(const DVector& v1, const DVector& v2) return t;
{
  t = v1;
  t -= v2;
}

DVector operator*(const DVector& v1, const DVector& v2) return t;
{
  t = v1;
  t *= v2;
}

DVector operator*(DVector& dv, DMatrix& dm) return t;
{
  if(dv.dim()!=dm.nrow())
    err_quit("operator*: vectors not same length %d (%d,%d)",
                dv.dim(),dm.nrow(),dm.ncol());

  t.resize(dm.ncol());

  for(int i=0; i < t.dim(); i++) {
    double tt=0;
    for(int j=0; j < dv.dim(); j++) tt += dv[j]*dm[j][i];
    t[i] = tt;
    }
}

DVector operator*(DMatrix& dm, DVector& dv) return t;
{
  if(dv.dim()!=dm.ncol())
    err_quit("operator*: vectors not same length %d (%d,%d)",
                dv.dim(),dm.nrow(),dm.ncol());

  t.resize(dm.nrow());

  for(int i=0; i < t.dim(); i++) {
    double tt=0;
    for (int j=0; j < dv.dim(); j++) tt += dm(i,j)*dv[j];
    t[i] = tt;
    }
}

DMatrix ccross(const DVector& v1, const DVector& v2)
  return t(v1.dim(),v2.dim());
{
  for(int i=0; i < v1.dim(); i++)
    for(int j=0; j < v2.dim(); j++)
      t(i,j) = v1(i)*v2(j);
}

#else /* !__GNUC__ */

DVector operator+(const DVector& dv, const double s)
{
  DVector t=dv;
  t += s;
  return t;
  }

DVector operator+(double s, const DVector& dv)
{
  DVector t=dv;
  t += s;
  return t;
  }

DVector operator*(const DVector& dv, const double s)
{
  DVector t=dv;
  t *= s;
  return t;
  }

DVector operator*(double s, const DVector& dv)
{
  DVector t=dv;
  t *= s;
  return t;
  }

DVector operator+(const DVector& v1, const DVector& v2)
{
  DVector t=v1;
  t += v2;
  return t;
  }

DVector operator-(const DVector& v1, const DVector& v2)
{
  DVector t=v1;
  t -= v2;
  return t;
  }

DVector operator*(const DVector& v1, const DVector& v2)
{
  DVector t=v1;
  t *= v2;
  return t;
  }

DVector operator*(DVector& dv, DMatrix& dm)
{
  if(dv.dim()!=dm.nrow())
    err_quit("operator*: vectors not same length %d (%d,%d)",
                dv.dim(),dm.nrow(),dm.ncol());

  DVector t(dm.ncol());

  for(int i=0; i < t.dim(); i++) {
    double tt=0;
    for(int j=0; j < dv.dim(); j++) tt += dv[j]*dm[j][i];
    t[i] = tt;
    }

  return t;
  }

DVector operator*(DMatrix& dm, DVector& dv)
{
  if(dv.dim()!=dm.ncol())
    err_quit("operator*: vectors not same length %d (%d,%d)",
                dv.dim(),dm.nrow(),dm.ncol());

  DVector t(dm.nrow());

  for(int i=0; i < t.dim(); i++) {
    double tt=0;
    for (int j=0; j < dv.dim(); j++) tt += dm(i,j)*dv[j];
    t[i] = tt;
    }

  return t;
  }

DMatrix ccross(const DVector& v1, const DVector& v2)
{
  DMatrix t(v1.dim(),v2.dim());

  for(int i=0; i < v1.dim(); i++)
    for(int j=0; j < v2.dim(); j++)
      t(i,j) = v1(i)*v2(j);

  return t;
}

#endif /* __GNUC__ */
