//
// transform.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <sstream>

#include <util/misc/formio.h>
#include <math/scmat/matrix.h>
#include <math/scmat/repl.h>
#include <chemistry/qc/basis/transform.h>

using namespace std;
using namespace sc;

#undef DEBUG

////////////////////////////////////////////////////////////////////////////
// Utility classes and routines to generate cartesian to pure transformation
// matrices.

class SafeUInt {
  private:
    unsigned long i_;
  public:
    SafeUInt(): i_(0) {}
    SafeUInt(unsigned long i): i_(i) {}
    void error() const {
      ExEnv::errn() << "SafeUInt: integer size exceeded" << endl;
      abort();
    }
    SafeUInt &operator ++ () { i_++; return *this; }
    SafeUInt &operator ++ (int) { i_++; return *this; }
    operator double() const { return i_; }
    operator unsigned long() const { return i_; }
    SafeUInt &operator =(const SafeUInt &i) { i_ = i.i_; return *this; }
    int operator > (const SafeUInt& i) const { return i_>i.i_; }
    int operator >= (const SafeUInt& i) const { return i_>=i.i_; }
    int operator < (const SafeUInt& i) const { return i_<i.i_; }
    int operator <= (const SafeUInt& i) const { return i_<=i.i_; }
    int operator == (const SafeUInt& i) const { return i_==i.i_; }
    int operator == (unsigned long i) const { return i_==i; }
    int operator != (const SafeUInt& i) const { return i_!=i.i_; }
    SafeUInt operator / (const SafeUInt& i) const { return SafeUInt(i_/i.i_); }
    SafeUInt operator % (const SafeUInt& i) const { return SafeUInt(i_%i.i_); }
    SafeUInt operator * (unsigned long i) const
    {
      unsigned long tmp = i_*i;
      if (tmp/i != i_ || tmp%i != 0) {
          error();
        }
      return SafeUInt(tmp);
    }
    SafeUInt operator * (const SafeUInt& i) const
    {
      return this->operator*(i.i_);
    }
    SafeUInt &operator = (unsigned long i) { i_ = i; return *this; }
    SafeUInt &operator *= (unsigned long i) { *this = *this*i; return *this; }
    SafeUInt &operator /= (unsigned long i) { i_ /= i; return *this; }
};

// there ordering here is arbitrary and doesn't have to match the
// basis set ordering
static inline int ncart(int l) { return (l>=0)?((((l)+2)*((l)+1))>>1):0; }
static inline int npure(int l) { return 2*l+1; }
static inline int icart(int a, int b, int c)
{
  return (((((a+b+c+1)<<1)-a)*(a+1))>>1)-b-1;
}

#define USE_OLD_SOLIDHARM_ORDERING 0
#if USE_OLD_SOLIDHARM_ORDERING
static inline int ipure(int l, int m) { return m<0?2*-m:(m==0?0:2*m-1); }
#else  // CCA ordering
static inline int ipure(int l, int m) { return l+m; }
#endif

static inline int local_abs(int i) { return i<0? -i:i; }

SafeUInt
binomial(int n, int c1)
{
  SafeUInt num = 1;
  SafeUInt den = 1;
  int c2 = n - c1;
  int i;
  for (i=c2+1; i<=n; i++) {
      num *= i;
    }
  for (i=2; i<=c1; i++) {
      den *= i;
    }
  return num/den;
}

SafeUInt
fact(int n)
{
  SafeUInt r = 1;
  for (int i=2; i<=n; i++) {
      r *= i;
    }
  return r;
}

// compute nnum!/nden!, nden <= nnum
SafeUInt
factoverfact(int nnum,int nden)
{
  SafeUInt r = 1;
  for (int i=nden+1; i<=nnum; i++) {
      r *= i;
    }
  return r;
}

SafeUInt
factfact(int n)
{
  SafeUInt result;
  int i;

  result = 1;
  if (n&1) {
      for (i=3; i<=n; i+=2) {
          result *= i;
        }
    }
  else {
      for (i=2; i<=n; i+=2) {
          result *= i;
        }
    }
  return result;
}

void
reduce(SafeUInt &num, SafeUInt &den)
{
  if (num > den) {
      for (SafeUInt i=2; i<=den;) {
          if (num%i == 0UL && den%i == 0UL) {
              num /= i;
              den /= i;
            }
          else i++;
        }
    }
  else {
      for (SafeUInt i=2; i<=num;) {
          if (num%i == 0UL && den%i == 0UL) {
              num /= i;
              den /= i;
            }
          else i++;
        }
    }
}

SafeUInt
powll(SafeUInt n, unsigned long p)
{
  SafeUInt result = 1;
  for (unsigned long i=0; i<p; i++) result *= n;
  return result;
}

////////////////////////////////////////////////////////////////////////////
// SphericalTransform class

SphericalTransform::SphericalTransform()
{
  n_ = 0;
  l_ = 0;
  subl_ = 0;
}

SphericalTransform::SphericalTransform(int l, int subl) : l_(l)
{
  n_ = 0;
  if (subl == -1) subl_ = l;
  else subl_ = subl;
}

static void
solidharmcontrib(int sign,
                 const SafeUInt &bin,const SafeUInt &den,
                 SafeUInt norm2num,SafeUInt norm2den,
                 int r2,int x,int y,int z,
                 const RefSCMatrix &coefmat, int pureindex)
{
  if (r2>0) {
      solidharmcontrib(sign,bin,den,norm2num,norm2den,r2-1,x+2,y,z,
                       coefmat,pureindex);
      solidharmcontrib(sign,bin,den,norm2num,norm2den,r2-1,x,y+2,z,
                       coefmat,pureindex);
      solidharmcontrib(sign,bin,den,norm2num,norm2den,r2-1,x,y,z+2,
                       coefmat,pureindex);
    }
  else {
      double coef = sign*double(bin)/double(den);
      double norm = sqrt(double(norm2num)/double(norm2den));
      coefmat->accumulate_element(icart(x,y,z), pureindex, coef*norm);
#ifdef DEBUG
      ExEnv::outn().form("    add(%d,%d,%d, % 4ld.0",
                x,y,z, sign*long(bin));
      if (den!=1) {
          ExEnv::outn().form("/%-4.1f",double(den));
        }
      else {
          ExEnv::outn().form("     ");
        }
      if (norm2num != 1 || norm2den != 1) {
          ExEnv::outn().form(" * sqrt(%ld.0/%ld.0)",
                    long(norm2num), long(norm2den));
        }
      ExEnv::outn().form(", i);");
      ExEnv::outn() << endl;
#endif
    }
}

// l is the total angular momentum
// m is the z component
// r2 is the number of factors of r^2 that are included
static void
solidharm(unsigned int l, int m, unsigned int r2, RefSCMatrix coefmat)
{
  int pureindex = ipure(l,m);
  for (unsigned int i=1; i<=r2; i++) pureindex += npure(l+2*i);
  
  unsigned int absm = local_abs(m);

  // the original norm2num and norm2den computation overflows 32bits for l=7
  //SafeUInt norm2num = factoverfact(l+absm,l-absm);
  //if (m != 0) norm2num *= 2;
  //SafeUInt normden = factfact(2*absm)*binomial(l,absm);
  //SafeUInt norm2den = normden*norm2den;
  //reduce(norm2num,norm2den);

  // this overflows 32bits for l=9
  SafeUInt norm2num = factoverfact(l+absm,l);
  SafeUInt norm2den = factoverfact(l,l-absm);
  reduce(norm2num,norm2den);
  norm2num *= fact(absm);
  norm2den *= factfact(2*absm);
  reduce(norm2num,norm2den);
  norm2num *= fact(absm);
  norm2den *= factfact(2*absm);
  if (m != 0) norm2num *= 2;
  reduce(norm2num,norm2den);

#ifdef DEBUG
  ExEnv::outn().form("    // l=%2d m=% 2d",l,m);
  ExEnv::outn() << endl;
#endif
  for (unsigned int t=0; t <= (l - absm)/2; t++) {
      for (unsigned int u=0; u<=t; u++) {
          int v2m;
          if (m >= 0) v2m = 0;
          else v2m = 1;
          for (unsigned int v2 = v2m; v2 <= absm; v2+=2) {
              int x = 2*t + absm - 2*u - v2;
              int y = 2*u + v2;
              int z = l - x - y;
              SafeUInt bin = binomial(l,t)
                               *binomial(l-t,absm+t)
                               *binomial(t,u)
                               *binomial(absm,v2);
              SafeUInt den = powll(4,t);
              int sign;
              if ((t + (v2-v2m)/2)%2) sign = -1;
              else sign = 1;
              reduce(bin,den);
              solidharmcontrib(sign,bin,den,norm2num,norm2den,
                               r2,x,y,z,coefmat,pureindex);
            }
        }
    }
#ifdef DEBUG
  ExEnv::outn() << "    i++;" << endl;
#endif
}

static void
solidharm(int l, const RefSCMatrix &coefmat)
{
  solidharm(l,0,0,coefmat);
  for (int m=1; m<=l; m++) {
      solidharm(l, m,0,coefmat);
      solidharm(l,-m,0,coefmat);
    }
  for (int r=2; r<=l; r+=2) {
      solidharm(l-r,0,r/2,coefmat);
      for (int m=1; m<=l-r; m++) {
          solidharm(l-r, m,r/2,coefmat);
          solidharm(l-r,-m,r/2,coefmat);
        }
    }

#ifdef DEBUG
  ExEnv::outn() << coefmat;
#endif
#if 0
  {
      ostringstream oss;
      oss << "cart->sph.harm coefficients l=" << l << endl;
      coefmat.print(oss.str().c_str());
  }
#endif
}

void
SphericalTransform::init()
{
  Ref<SCMatrixKit> matrixkit = new ReplSCMatrixKit;
  RefSCDimension cartdim(new SCDimension(ncart(l_)));
  RefSCMatrix coefmat(cartdim,cartdim,matrixkit);
  coefmat->assign(0.0);

  solidharm(l_,coefmat);

#ifdef DEBUG
  ExEnv::outn() << scprintf("---> generating l=%d subl=%d", l_, subl_) << endl;
#endif

  int pureoffset = 0;
  for (int i=1; i<=(l_-subl_)/2; i++) pureoffset += npure(subl_+2*i);

  for (int p=0; p<npure(subl_); p++) {
    for (int a=0; a<=l_; a++) {
      for (int b=0; (a+b)<=l_; b++) {
        int c = l_ - a - b;
        int cart = icart(a,b,c);
        double coef = coefmat->get_element(cart,p+pureoffset);
        if (fabs(coef) > DBL_EPSILON) {
          add(a,b,c, coef, p);
#ifdef DEBUG
          ExEnv::outn() << scprintf("---> add(%d,%d,%d, %12.8f, %d)",
                           a,b,c,coef,p) << endl;
#endif
        }
      }
    }
  }
}

SphericalTransform::~SphericalTransform()
{
  for (auto* c: components_) {
    delete c;
  }
}

void
SphericalTransform::add(int a, int b, int c, double coef, int pureindex)
{
  int i;

  n_++;
  components_.resize(n_);
  components_.back() = new_component();
  components_.back()->init(a, b, c, coef, pureindex);
}

///////////////////////////////////////////////////////////////////////////

ISphericalTransform::ISphericalTransform() :
  SphericalTransform()
{
}

ISphericalTransform::ISphericalTransform(int l,int subl) :
  SphericalTransform(l,subl)
{
}

void
ISphericalTransform::init()
{
  Ref<SCMatrixKit> matrixkit = new ReplSCMatrixKit;
  RefSCDimension cartdim(new SCDimension(ncart(l_)));
  RefSCMatrix coefmat(cartdim,cartdim,matrixkit);
  coefmat->assign(0.0);

  solidharm(l_,coefmat);

  coefmat->invert_this();
  coefmat->transpose_this();

#ifdef DEBUG
  ExEnv::outn() << scprintf("---> IST: generating l=%d subl=%d", l_, subl_) << endl;
#endif

  int pureoffset = 0;
  for (int i=1; i<=(l_-subl_)/2; i++) pureoffset += npure(subl_+2*i);

  for (int p=0; p<npure(subl_); p++) {
    for (int a=0; a<=l_; a++) {
      for (int b=0; (a+b)<=l_; b++) {
        int c = l_ - a - b;
        int cart = icart(a,b,c);
        double coef = coefmat->get_element(cart,p+pureoffset);
        if (fabs(coef) > DBL_EPSILON) {
          add(a,b,c, coef, p);
#ifdef DEBUG
          ExEnv::outn() << scprintf("---> IST: add(%d,%d,%d, %12.8f, %d)",
                           a,b,c,coef,p) << endl;
#endif
        }
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////

SphericalTransformIter::SphericalTransformIter()
{
  transform_=0;
}

SphericalTransformIter::SphericalTransformIter(const SphericalTransform*t)
{
  transform_ = t;
}

////////////////////////////////////////////////////////////////////////////
// SphericalTransformComponent class

SphericalTransformComponent::~SphericalTransformComponent() {}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
