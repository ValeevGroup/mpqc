//
// bend.cc
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

/* bend.cc -- implementation of the bending simple internal coordinate class
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

static ClassDesc BendSimpleCo_cd(
  typeid(BendSimpleCo),"BendSimpleCo",1,"public SimpleCo",
  create<BendSimpleCo>, create<BendSimpleCo>, create<BendSimpleCo>);
SimpleCo_IMPL(BendSimpleCo)

BendSimpleCo::BendSimpleCo() : SimpleCo(3) {}

BendSimpleCo::BendSimpleCo(const BendSimpleCo& s)
  : SimpleCo(3)
{
  *this=s;
}

BendSimpleCo::BendSimpleCo(const char *refr, int a1, int a2, int a3)
  : SimpleCo(3,refr)
{
  atoms[0]=a1; atoms[1]=a2; atoms[2]=a3;
}

BendSimpleCo::BendSimpleCo(const Ref<KeyVal> &kv)
  : SimpleCo(kv,3)
{
}

BendSimpleCo::~BendSimpleCo()
{
}

BendSimpleCo&
BendSimpleCo::operator=(const BendSimpleCo& s)
{
  if(label_) delete[] label_;
  label_=new char[strlen(s.label_)+1]; strcpy(label_,s.label_);
  atoms[0]=s.atoms[0]; atoms[1]=s.atoms[1]; atoms[2]=s.atoms[2];
  return *this;
}

double
BendSimpleCo::calc_intco(Molecule& m, double *bmat, double coeff)
{
  SCVector3 u1, u2;
  int a=atoms[0]-1; int b=atoms[1]-1; int c=atoms[2]-1;

  SCVector3 ra(m.r(a));
  SCVector3 rb(m.r(b));
  SCVector3 rc(m.r(c));

  u1 = ra-rb;
  u1.normalize();
  u2 = rc-rb;
  u2.normalize();

  double co=u1.dot(u2);

  value_=acos(co);

  if(bmat) {
    double uu,ww,vv;
    double si=s2(co);
    double r1i, r2i;
    if (si > 1.0e-4) {
        r1i = 1.0/(si*ra.dist(rb));
        r2i = 1.0/(si*rc.dist(rb));
      }
    else {r1i = 0.0; r2i = 0.0;}
#if OLD_BMAT
    r1i /= bohr;
    r2i /= bohr;
#endif    
    for (int j=0; j < 3; j++) {
      uu = (co*u1[j]-u2[j])*r1i;
      ww = (co*u2[j]-u1[j])*r2i;
      vv = -uu-ww;
      bmat[a*3+j] += coeff*uu;
      bmat[b*3+j] += coeff*vv;
      bmat[c*3+j] += coeff*ww;
    }
  }

  return value_;
}

double
BendSimpleCo::calc_force_con(Molecule& m)
{
  int a=atoms[1]-1; int b=atoms[0]-1; int c=atoms[2]-1;

  double rad_ab =   m.atominfo()->atomic_radius(m.Z(a))
                  + m.atominfo()->atomic_radius(m.Z(b));

  double rad_ac =   m.atominfo()->atomic_radius(m.Z(a))
                  + m.atominfo()->atomic_radius(m.Z(c));

  SCVector3 ra(m.r(a));
  SCVector3 rb(m.r(b));
  SCVector3 rc(m.r(c));

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
BendSimpleCo::ctype() const
{
  return "BEND";
}

double
BendSimpleCo::radians() const
{
  return value_;
}

double
BendSimpleCo::degrees() const
{
  return value_*rtd;
}

double
BendSimpleCo::preferred_value() const
{
  return value_*rtd;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
