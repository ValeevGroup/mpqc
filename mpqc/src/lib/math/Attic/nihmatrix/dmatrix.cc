
#include "nihmatrix.h"
#include "lmath.h"
#include <iostream.h>
#include <util/keyval/keyval.h>

DescribedClass_REF_def(DMatrix);

#define CLASSNAME DMatrix
#define PARENTS virtual public SavableState
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
DMatrix::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}

DMatrix::DMatrix(int s1, int s2) : n1(0), n2(0), d(0), dp(0)
{
  resize(s1,s2);
  }

DMatrix::DMatrix(int s1, int s2, double *dnew) : n1(s1), n2(s2), d(dnew), dp(0)
{
  dp = new double*[s1];
  for(int i=0; i < s1; i++) dp[i]=&d[i*s2];
  }

DMatrix::DMatrix() : n1(0), n2(0), d(0), dp(0) {}
DMatrix::DMatrix(const DMatrix& da) : n1(0), n2(0), d(0), dp(0) { *this=da; }

DMatrix::~DMatrix() { init(); }

void DMatrix::init()
{
  if(d) delete[] d; d=0;
  if(dp) delete[] dp; dp=0;
  n1=n2=0;
  }

void DMatrix::resize(int s1, int s2)
{
  if(!s1*s2) {
    init();
    return;
    }

  if(s1*s2 != n1*n2) {
    if(d) delete[] d; d = new double[s1*s2];
    }

  if(s1!=n1 || s2!=n2) {
    if(dp) delete[] dp; dp = new double*[s1];
    for(int i=0; i < s1; i++) dp[i]=&d[i*s2];
    }

  n1=s1; n2=s2;
  }

DMatrix& DMatrix::operator=(const DMatrix& dm)
{
  resize(dm.n1,dm.n2);

  for(int i=0; i < n1*n2; i++) d[i]=dm.d[i];
  return *this;
  }

DMatrix& DMatrix::operator+=(const double s)
{
  for(int i=0; i < n1*n2; i++) d[i] += s;
  return *this;
  }
      
DMatrix& DMatrix::operator*=(const double s)
{
  for(int i=0; i < n1*n2; i++) d[i] *= s;
  return *this;
  }
      
DMatrix& DMatrix::operator+=(const DMatrix& dm)
{
  if(n1!=dm.n1 || n2!=dm.n2)
    err_quit("DMatrix::operator+=: matrices are wrong size");

  for(int i=0; i < n1*n2; i++) d[i] += dm.d[i];
  return *this;
  }

DMatrix& DMatrix::operator-=(const DMatrix& dm)
{
  if(n1!=dm.n1 || n2!=dm.n2)
    err_quit("DMatrix::operator-=: matrices are wrong size");

  for(int i=0; i < n1*n2; i++) d[i] -= dm.d[i];
  return *this;
  }

DMatrix& DMatrix::operator*=(const DMatrix& dm)
{
  if(n1!=dm.n1 || n2!=dm.n2)
    err_quit("DMatrix::operator*=: matrices are wrong size");

  for(int i=0; i < n1*n2; i++) d[i] *= dm.d[i];
  return *this;
  }

///////////////////////////////////////////////////////////////

#ifdef __GNUC__
DMatrix DMatrix::transpose() const return t(n2,n1);
{
  for(int i=0; i < n1; i++)
    for(int j=0; j < n2; j++) t[j][i]=dp[i][j];
}

DMatrix DMatrix::inverse() const return t;
{
  t = *this;
  t.invert();
}
#else
DMatrix DMatrix::transpose() const
{
  DMatrix t(n2,n1);
  for(int i=0; i < n1; i++)
    for(int j=0; j < n2; j++) t[j][i]=dp[i][j];
  return t;
  }

DMatrix DMatrix::inverse() const
{
  DMatrix t= *this;
  t.invert();
  return t;
  }
#endif

void DMatrix::transpose_IP()
{
  double *dx=new double[n1*n2];
  for(int i=0; i < n1*n2; i++) dx[i]=d[i];
  resize(n2,n1);
  for(i=0; i < n1; i++)
    for(int j=0; j < n2; j++) dp[i][j] = dx[j*n1+i];
  delete[] dx;
}

double DMatrix::invert()
{
  return mathqc_invert(this);
}

double DMatrix::solve_lin(DVector& b)
{
  return mathqc_solve_lin(*this,b);
}

#ifdef __GNUC__
DMatrix DMatrix::eigenvectors() return t(n1,n2);
{
  DVector e(n1);

  mathqc_diag(*this,e,t,1);
  return t;
}

DVector DMatrix::eigenvalues() return e(n1);
{
  DMatrix t(n1,n2);

  mathqc_diag(*this,e,t,0);
  return e;
}
#else
DMatrix DMatrix::eigenvectors()
{
  DMatrix t(n1,n2);
  DVector e(n1);

  mathqc_diag(*this,e,t,1);
  return t;
}

DVector DMatrix::eigenvalues()
{
  DMatrix t(n1,n2);
  DVector e(n1);

  mathqc_diag(*this,e,t,0);
  return e;
}
#endif

void DMatrix::diagonalize(DVector& evals, DMatrix& evecs, double tol)
{
  mathqc_diag(*this,evals,evecs,1,tol);
}

///////////////////////////////////////////////////////////////

DMatrix::DMatrix(KeyVal&k):
  n1(0),
  n2(0),
  d(0),
  dp(0)
{
  int s1,s2;
  s1 = k.count();
  if (k.error() != KeyVal::OK) {
      s1 = s2 = 0;
    }
  else {
      s2 = k.count("0");
    }
  resize(s1,s2);
  for (int i=0; i<n1; i++) {
      for (int j=0; j<n2; j++) {
          dp[i][j] = k.doublevalue(i,j);
        }
    }
}

void DMatrix::save_data_state(StateOut& so)
{
  so.put(n1); so.put(n2); so.put(d,n1*n2);
}

DMatrix::DMatrix(StateIn& si):
  SavableState(si,class_desc_)
{
  si.get(n1); si.get(n2); si.get(d);
  int refnum;
  dp = new double*[n1];
  for(int i=0; i < n1; i++) dp[i]=&d[i*n2];
}

const double DMatrix::maxval() const
{
  double t=0,tt;
  for(int i=0; i < n1*n2; i++) if((tt=fabs(d[i]))>t) t=tt;
  return t;
  }

void DMatrix::print(ostream& os) const { print(0,os); }

void DMatrix::print(const char *title, ostream& os, int prec) const
{
  int ii,jj,kk,nn,ll;
  int i,j;
  int lwidth,width;
  double max=this->maxval();

  max=(max==0.0)?1.0:log10(max);
  if(max < 0.0) max=1.0;

  lwidth = prec+5+(int) max; width = 75/lwidth;

  os.setf(ios::fixed,ios::floatfield); os.precision(prec);
  os.setf(ios::right,ios::adjustfield);

  if(title) os << "\n" << title << "\n";
  else os << "\n";

  if(n1==0 || n2==0) { os << " empty matrix\n"; return; }

  for(ii=jj=0;;) {
    ii++; jj++; kk=width*jj;
    nn=(n2>kk)?kk:n2;
    ll=2*(nn-ii+1)+1;

 // print column indices
    for(i=ii; i <= nn; i++) { os.width(lwidth); os << i; }
    os << "\n";

 // print the rows
    for(i=0; i < n1 ; i++) {
      os.width(5); os << i+1;
      for(j=ii-1; j < nn; j++) { os.width(lwidth); os << dp[i][j]; }
      os << "\n";
      }
    os << "\n";

    if(n2<=kk) { os.flush(); return; }
    ii=kk;
    }
  }

////////////////////////////////////////////////////////////////////

#ifdef __GNUC__
DMatrix operator+(DMatrix& m1, DMatrix& m2) return t;
{
  if(m1.ncol()!=m2.ncol() || m1.nrow() != m2.nrow())
    err_quit("operator+: matrices are wrong size, (%d,%d) (%d,%d)",
              m1.nrow(),m1.ncol(),m2.nrow(),m2.ncol());

  t = m1;
  t += m2;
}

DMatrix operator-(DMatrix& m1, DMatrix& m2) return t;
{
  if(m1.ncol()!=m2.ncol() || m1.nrow() != m2.nrow())
    err_quit("operator-: matrices are wrong size, (%d,%d) (%d,%d)",
              m1.nrow(),m1.ncol(),m2.nrow(),m2.ncol());

  t = m1;
  t -= m2;
}

DMatrix operator*(DMatrix& m1, DMatrix& m2) return t;
{
  if(m1.ncol()!=m2.nrow())
    err_quit("operator*: matrices are wrong size, (%d,%d) (%d,%d)",
              m1.nrow(),m1.ncol(),m2.nrow(),m2.ncol());

  t.resize(m1.nrow(),m2.ncol());
  mathqc_mxm(m1,0,m2,0,t,0,m1.nrow(),m1.ncol(),m2.ncol(),0);
}

DMatrix operator*(DMatrix& m, double s) return t;
{
  t = m;
  t *= s;
}

DMatrix operator*(double s, DMatrix& m) return t;
{
  t = m;
  t *= s;
}

#else /* !__GNUC__ */

DMatrix operator+(DMatrix& m1, DMatrix& m2)
{
  if(m1.ncol()!=m2.ncol() || m1.nrow() != m2.nrow())
    err_quit("operator+: matrices are wrong size, (%d,%d) (%d,%d)",
              m1.nrow(),m1.ncol(),m2.nrow(),m2.ncol());

  DMatrix t=m1;
  t += m2;
  return t;
  }

DMatrix operator-(DMatrix& m1, DMatrix& m2)
{
  if(m1.ncol()!=m2.ncol() || m1.nrow() != m2.nrow())
    err_quit("operator-: matrices are wrong size, (%d,%d) (%d,%d)",
              m1.nrow(),m1.ncol(),m2.nrow(),m2.ncol());

  DMatrix t=m1;
  t -= m2;
  return t;
  }

DMatrix operator*(DMatrix& m1, DMatrix& m2)
{
  if(m1.ncol()!=m2.nrow())
    err_quit("operator*: matrices are wrong size, (%d,%d) (%d,%d)",
              m1.nrow(),m1.ncol(),m2.nrow(),m2.ncol());

  DMatrix t(m1.nrow(),m2.ncol());
  mathqc_mxm(m1,0,m2,0,t,0,m1.nrow(),m1.ncol(),m2.ncol(),0);
  return t;
  }

DMatrix operator*(DMatrix& m, double s)
{
  DMatrix t=m;
  t *= s;
  return t;
  }

DMatrix operator*(double s, DMatrix& m)
{
  DMatrix t=m;
  t *= s;
  return t;
  }

#endif
