//
// tors.cc
//
// Modifications are
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
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

/* tors.cc -- implementation of the torsion internal coordinate class
 *
 *      THIS SOFTWARE FITS THE DESCRIPTION IN THE U.S. COPYRIGHT ACT OF A
 *      "UNITED STATES GOVERNMENT WORK".  IT WAS WRITTEN AS A PART OF THE
 *      AUTHOR'S OFFICIAL DUTIES AS A GOVERNMENT EMPLOYEE.  THIS MEANS IT
 *      CANNOT BE COPYRIGHTED.  THIS SOFTWARE IS FREELY AVAILABLE TO THE
 *      PUBLIC FOR USE WITHOUT A COPYRIGHT NOTICE, AND THERE ARE NO
 *      RESTRICTIONS ON ITS USE, NOW OR SUBSEQUENTLY.
 *
 *  Author:
 *      E. T. Seidl
 *      Bldg. 12A, Rm. 2033
 *      Computer Systems Laboratory
 *      Division of Computer Research and Technology
 *      National Institutes of Health
 *      Bethesda, Maryland 20892
 *      Internet: seidl@alw.nih.gov
 *      February, 1993
 */

#include <string.h>
#include <math.h>

#include <chemistry/molecule/simple.h>
#include <chemistry/molecule/localdef.h>

using namespace sc;

static ClassDesc TorsSimpleCo_cd(
  typeid(TorsSimpleCo),"TorsSimpleCo",1,"public SimpleCo",
  create<TorsSimpleCo>, create<TorsSimpleCo>, create<TorsSimpleCo>);
SimpleCo_IMPL(TorsSimpleCo)


TorsSimpleCo::TorsSimpleCo() : SimpleCo(4) {}

TorsSimpleCo::TorsSimpleCo(const TorsSimpleCo& s)
  : SimpleCo(4)
{
  *this=s;
}

TorsSimpleCo::TorsSimpleCo(const char *refr, int a1, int a2, int a3, int a4)
  : SimpleCo(4,refr)
{
  atoms[0]=a1; atoms[1]=a2; atoms[2]=a3; atoms[3]=a4;
}

TorsSimpleCo::~TorsSimpleCo()
{
}

TorsSimpleCo::TorsSimpleCo(const Ref<KeyVal> &kv):
  SimpleCo(kv,4)
{
}

TorsSimpleCo&
TorsSimpleCo::operator=(const TorsSimpleCo& s)
{
  if(label_) delete[] label_;
  label_=new char[strlen(s.label_)+1];
  strcpy(label_,s.label_);
  atoms[0]=s.atoms[0]; atoms[1]=s.atoms[1]; atoms[2]=s.atoms[2];
  atoms[3]=s.atoms[3];
  return *this;
}

double
TorsSimpleCo::calc_intco(Molecule& m, double *bmat, double coeff)
{
  int a=atoms[0]-1; int b=atoms[1]-1; int c=atoms[2]-1; int d=atoms[3]-1;
  SCVector3 u1,u2,u3,z1,z2;

  SCVector3 ra(m.r(a));
  SCVector3 rb(m.r(b));
  SCVector3 rc(m.r(c));
  SCVector3 rd(m.r(d));

  u1 = ra-rb;
  u1.normalize();
  u2 = rc-rb;
  u2.normalize();
  u3 = rc-rd;
  u3.normalize();

  z1 = u1.perp_unit(u2);
  z2 = u3.perp_unit(u2);

  double co=z1.dot(z2);
  u1[0]=z1[1]*z2[2]-z1[2]*z2[1];
  u1[1]=z1[2]*z2[0]-z1[0]*z2[2];
  u1[2]=z1[0]*z2[1]-z1[1]*z2[0];
  double co2=u1.dot(u2);

  if (co < -1.0) co= -1.0;
  if (co > 1.0) co = 1.0;

  // save the old value of the torsion so we can make sure the discontinuity
  // at -pi/2 doesn't bite us

  double oldval = -value_;  

  value_=(co2<0) ? -acos(-co) : acos(-co);

  // ok, we want omega between 3*pi/2 and -pi/2, so if omega is > pi/2
  // (omega is eventually -omega), then knock 2pi off of it
  if(value_ > pih) value_ -= tpi;

  // the following tests to see if the new coordinate has crossed the
  // 3pi/2 <--> -pi/2 boundary...if so, then we add or subtract 2pi as
  // needed to prevent the transformation from internals to cartesians
  // from blowing up
  while(oldval-value_ > (pi + 1.0e-8)) value_ += tpi;
  while(oldval-value_ < -(pi + 1.0e-8)) value_ -= tpi;

  value_ = -value_;

  if (bmat) {
    double uu,vv,ww,zz;
    u1 = ra-rb;
    u1.normalize();
    u2 = rc-rb;
    u2.normalize();
    u3 = rc-rd;
    u3.normalize();
    z1 = u1.perp_unit(u2);
    z2 = u3.perp_unit(u2);
    co=u1.dot(u2); double si=s2(co);
    co2=u2.dot(u3); double si2=s2(co2);
    double r1 = ra.dist(rb);
    double r2 = rc.dist(rb);
    double r3 = rc.dist(rd);
#if OLD_BMAT
    r1 *= bohr;
    r2 *= bohr;
    r3 *= bohr;
#endif    
    for (int j=0; j < 3; j++) {
      if (si > 1.0e-5) uu = z1[j]/(r1*si);
      else uu = 0.0;
      if (si2 > 1.0e-5) zz = z2[j]/(r3*si2);
      else zz = 0.0;
      vv = (r1*co/r2-1.0)*uu-zz*r3*co2/r2;
      ww = -uu-vv-zz;
      bmat[a*3+j] += coeff*uu;
      bmat[b*3+j] += coeff*vv;
      bmat[c*3+j] += coeff*ww;
      bmat[d*3+j] += coeff*zz;
    }
  }

  return value_;
}

double
TorsSimpleCo::calc_force_con(Molecule& m)
{
  int a=atoms[1]-1; int b=atoms[2]-1;

  double rad_ab =   m.atominfo()->atomic_radius(m.Z(a))
                  + m.atominfo()->atomic_radius(m.Z(b));

  SCVector3 ra(m.r(a));
  SCVector3 rb(m.r(b));

  double r_ab = ra.dist(rb);

  double k = 0.0015 + 14.0*pow(1.0,0.57)/pow((rad_ab*r_ab),4.0) *
                           exp(-2.85*(r_ab-rad_ab));

#if OLD_BMAT  
  // return force constant in mdyn*ang/rad^2
  return k*4.359813653;
#else
  return k;
#endif  
}

const char *
TorsSimpleCo::ctype() const
{
  return "TORS";
}

double
TorsSimpleCo::radians() const
{
  return value_;
}

double
TorsSimpleCo::degrees() const
{
  return value_*rtd;
}

double
TorsSimpleCo::preferred_value() const
{
  return value_*rtd;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
