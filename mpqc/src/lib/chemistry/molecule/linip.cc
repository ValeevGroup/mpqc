//
// linip.cc
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

/* lin.cc -- implementation of the linear bending internal coordinate classes
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

static ClassDesc LinIPSimpleCo_cd(
  typeid(LinIPSimpleCo),"LinIPSimpleCo",1,"public SimpleCo",
  create<LinIPSimpleCo>, create<LinIPSimpleCo>, create<LinIPSimpleCo>);
SimpleCo_IMPL(LinIPSimpleCo)

LinIPSimpleCo::LinIPSimpleCo() : SimpleCo(3)
{
  u2[0] = 1.0; u2[1] = 0.0; u2[2] = 0.0;
}

LinIPSimpleCo::LinIPSimpleCo(const LinIPSimpleCo& s)
  : SimpleCo(3)
{
  *this=s;
}

LinIPSimpleCo::LinIPSimpleCo(const char *refr, int a1, int a2, int a3,
                             const SCVector3 &u)
  : SimpleCo(3,refr), u2(u)
{
  atoms[0]=a1; atoms[1]=a2; atoms[2]=a3;
  u2.normalize();
}

LinIPSimpleCo::~LinIPSimpleCo()
{
}

LinIPSimpleCo::LinIPSimpleCo(const Ref<KeyVal> &kv) :
  SimpleCo(kv,3)
{
  for (int i=0; i<3; i++) u2[i] = kv->doublevalue("u",i);
  u2.normalize();
}

LinIPSimpleCo&
LinIPSimpleCo::operator=(const LinIPSimpleCo& s)
{
  if(label_) delete[] label_;
  label_=new char[strlen(s.label_)+1];
  strcpy(label_,s.label_);
  atoms[0]=s.atoms[0]; atoms[1]=s.atoms[1]; atoms[2]=s.atoms[2];
  u2 = s.u2;
  return *this;
}

double
LinIPSimpleCo::calc_intco(Molecule& m, double *bmat, double coeff)
{
  int a=atoms[0]-1; int b=atoms[1]-1; int c=atoms[2]-1;
  SCVector3 u1,u3;

  SCVector3 ra(m.r(a));
  SCVector3 rb(m.r(b));
  SCVector3 rc(m.r(c));

  u1=ra-rb;
  u1.normalize();
  u3=rc-rb;
  u3.normalize();

  double co=u1.dot(u2);
  double co2=u3.dot(u2);

  value_ = pi-acos(co)-acos(co2);

  if (bmat) {
    double uu,ww,vv;
    SCVector3 z1,z2;
    z2 = u2.perp_unit(u1);
    z1 = u1.perp_unit(z2);
    z2 = u3.perp_unit(u2);
    u1 = z2.perp_unit(u3);
    double r1 = ra.dist(rb);
    double r2 = rc.dist(rb);
#if OLD_BMAT
    r1 *= bohr;
    r2 *= bohr;
#endif    
    for (int j=0; j < 3; j++) {
      uu=z1[j]/r1;
      ww=u1[j]/r2;
      vv = -uu-ww;
      bmat[a*3+j] += coeff*uu;
      bmat[b*3+j] += coeff*vv;
      bmat[c*3+j] += coeff*ww;
    }
  }

  return value_;
}

double
LinIPSimpleCo::calc_force_con(Molecule&m)
{
  int a=atoms[1]-1; int b=atoms[0]-1; int c=atoms[2]-1;

  SCVector3 ra(m.r(a));
  SCVector3 rb(m.r(b));
  SCVector3 rc(m.r(c));

  double rad_ab =   m.atominfo()->atomic_radius(m.Z(a))
                  + m.atominfo()->atomic_radius(m.Z(b));

  double rad_ac =   m.atominfo()->atomic_radius(m.Z(a))
                  + m.atominfo()->atomic_radius(m.Z(c));

  double r_ab = ra.dist(rb);
  double r_ac = ra.dist(rc);

  double k = 0.089 + 0.11/pow((rad_ab*rad_ac),-0.42) *
                           exp(-0.44*(r_ab+r_ac-rad_ab-rad_ac));

#if OLD_BMAT
  // return force constant in mdyn*ang/rad^2
  return k*4.359813653;
#else  
  return k;
#endif  
}

const char *
LinIPSimpleCo::ctype() const
{
  return "LINIP";
}

double
LinIPSimpleCo::radians() const
{
  return value_;
}

double
LinIPSimpleCo::degrees() const
{
  return value_*rtd;
}

double
LinIPSimpleCo::preferred_value() const
{
  return value_*rtd;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
