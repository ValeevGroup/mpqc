//
// shellrot.cc
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

#include <util/misc/formio.h>

#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/shellrot.h>
#include <chemistry/qc/basis/cartiter.h>
#include <chemistry/qc/basis/transform.h>

using namespace std;
using namespace sc;

void
ShellRotation::done() {
  if (r) {
    for (int i=0; i < n_; i++) {
      if (r[i]) delete[] r[i];
    }
    delete[] r;
    r=0;
  }
  n_=0;
}

ShellRotation::ShellRotation(int n) :
  n_(n),
  am_(0),
  r(0)
{
  if (n_) {
    r = new double*[n_];
    for (int i=0; i < n_; i++)
      r[i] = new double[n_];
  }
}

ShellRotation::ShellRotation(const ShellRotation& rot) :
  n_(0),
  am_(0),
  r(0)
{
  *this = rot;
}

ShellRotation::ShellRotation(int a, SymmetryOperation& so,
                             const Ref<Integral>& ints,
                             int pure) :
  n_(0),
  am_(0),
  r(0)
{
  if (a > 0 && pure)
    init_pure(a,so,ints);
  else
    init(a,so,ints);
}

ShellRotation::~ShellRotation()
{
  done();
}

ShellRotation&
ShellRotation::operator=(const ShellRotation& rot)
{
  done();

  n_ = rot.n_;
  am_ = rot.am_;

  if (n_ && rot.r) {
    r = new double*[n_];
    for (int i=0; i < n_; i++) {
      r[i] = new double[n_];
      memcpy(r[i],rot.r[i],sizeof(double)*n_);
    }
  }

  return *this;
}

// Compute the transformation matrices for general cartesian shells
// using the P (xyz) transformation matrix.  This is done as a
// matrix outer product, keeping only the unique terms.
// Written by clj...blame him
void
ShellRotation::init(int a, SymmetryOperation& so, const Ref<Integral>& ints)
{
  done();

  am_=a;

  if (a == 0) {
    n_ = 1;
    r = new double*[1];
    r[0] = new double[1];
    r[0][0] = 1.0;
    return;
  }
  
  CartesianIter *ip = ints->new_cartesian_iter(am_);
  RedundantCartesianIter *jp = ints->new_redundant_cartesian_iter(am_);
  
  CartesianIter& I = *ip;
  RedundantCartesianIter& J = *jp;
  int lI[3];
  int k, iI;
  
  n_ = I.n();
  r = new double*[n_];

  for (I.start(); I; I.next()) {
    r[I.bfn()] = new double[n_];
    memset(r[I.bfn()],0,sizeof(double)*n_);

    for (J.start(); J; J.next()) {
      double tmp = 1.0;

      for (k=0; k < 3; k++) {
        lI[k] = I.l(k);
      }
      
      for (k=0; k < am_; k++) {
        for (iI=0; lI[iI]==0; iI++);
        lI[iI]--;
        double contrib = so(J.axis(k),iI);
        tmp *= contrib;
      }

      r[I.bfn()][J.bfn()] += tmp;
    }
  }

  delete ip;
  delete jp;
}

// Compute the transformation matrices for general pure am
// by summing contributions from the cartesian components
// using the P (xyz) transformation matrix.  This is done as a
// matrix outer product, keeping only the unique terms.
void
ShellRotation::init_pure(int a, SymmetryOperation&so, const Ref<Integral>& ints)
{
  if (a < 1) {
    init(a,so,ints);
    return;
  }

  done();

  am_=a;
  
  SphericalTransformIter *ip = ints->new_spherical_transform_iter(am_);
  SphericalTransformIter *jp = ints->new_spherical_transform_iter(am_, 1);
  RedundantCartesianSubIter *kp = ints->new_redundant_cartesian_sub_iter(am_);
  
  SphericalTransformIter& I = *ip;
  SphericalTransformIter& J = *jp;
  RedundantCartesianSubIter& K = *kp;
  int lI[3];
  int m, iI;
  
  n_ = I.n();

  r = new double*[n_];
  for (m=0; m<n_; m++) {
      r[m] = new double[n_];
      memset(r[m],0,sizeof(double)*n_);
    }

  for (I.start(); I; I.next()) {
      for (J.start(); J; J.next()) {
          double coef = I.coef()*J.coef();
          double tmp = 0.0;
          for (K.start(J.a(), J.b(), J.c()); K; K.next()) {
              //printf("T(%d,%d) += %6.4f", I.bfn(), J.bfn(), coef);
              double tmp2 = coef;
              for (m=0; m < 3; m++) {
                  lI[m] = I.l(m);
                }
      
              for (m=0; m < am_; m++) {
                  for (iI=0; lI[iI]==0; iI++);
                  lI[iI]--;
                  //tmp2 *= so(iI,K.axis(m));
                  tmp2 *= so(K.axis(m),iI);
                  //printf(" * so(%d,%d) [=%4.2f]",
                  //       iI,K.axis(m),so(iI,K.axis(m)));
                }
              //printf(" = %8.6f\n", tmp2);
              tmp += tmp2;
            }
          r[I.bfn()][J.bfn()] += tmp;
        }
    }

  delete ip;
  delete jp;
  delete kp;
  
}

// returns the result of rot*this
ShellRotation
ShellRotation::operate(const ShellRotation& rot) const
{
  if (n_ != rot.n_) {
    ExEnv::err0() << indent
         << "ShellRotation::operate(): dimensions don't match" << endl
         << indent << scprintf("  %d != %d\n",rot.n_,n_);
    abort();
  }
  
  ShellRotation ret(n_);
  ret.am_ = am_;
  
  for (int i=0; i < n_; i++) {
    for (int j=0; j < n_; j++) {
      double t=0;
      for (int k=0; k < n_; k++)
        t += rot.r[i][k] * r[k][j];
      ret.r[i][j] = t;
    }
  }

  return ret;
}

ShellRotation
ShellRotation::transform(const ShellRotation& rot) const
{
  int i,j,k;

  if (rot.n_ != n_) {
    ExEnv::err0() << indent
         << "ShellRotation::transform(): dimensions don't match" << endl
         << indent << scprintf("%d != %d\n",rot.n_,n_);
    abort();
  }
  
  ShellRotation ret(n_), foo(n_);
  ret.am_ = foo.am_ = am_;

  // foo = r * d
  for (i=0; i < n_; i++) {
    for (j=0; j < n_; j++) {
      double t=0;
      for (k=0; k < n_; k++)
        t += rot.r[i][k] * r[k][j];
      foo.r[i][j] = t;
    }
  }

  // ret = (r*d)*r~ = foo*r~
  for (i=0; i < n_; i++) {
    for (j=0; j < n_; j++) {
      double t=0;
      for (k=0; k < n_; k++)
        t += foo.r[i][k]*rot.r[j][k];
      ret.r[i][j]=t;
    }
  }

  return ret;
}
    
double
ShellRotation::trace() const {
  double t=0;
  for (int i=0; i < n_; i++)
    t += r[i][i];
  return t;
}

void
ShellRotation::print() const
{
  for (int i=0; i < n_; i++) {
    ExEnv::out0() << indent << scprintf("%5d ",i+1);
    for (int j=0; j < n_; j++) {
      ExEnv::out0() << scprintf(" %10.7f",r[i][j]);
    }
    ExEnv::out0() << endl;
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
