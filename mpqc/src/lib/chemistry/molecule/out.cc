//
// out.cc
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

/* out.cc -- implementation of the out-of-plane internal coordinate class
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

static ClassDesc OutSimpleCo_cd(
  typeid(OutSimpleCo),"OutSimpleCo",1,"public SimpleCo",
  create<OutSimpleCo>, create<OutSimpleCo>, create<OutSimpleCo>);
SimpleCo_IMPL(OutSimpleCo)

OutSimpleCo::OutSimpleCo() : SimpleCo(4) {}

OutSimpleCo::OutSimpleCo(const OutSimpleCo& s)
  : SimpleCo(4)
{
  *this=s;
}

OutSimpleCo::OutSimpleCo(const char *refr, int a1, int a2, int a3, int a4)
  : SimpleCo(4,refr)
{
  atoms[0]=a1; atoms[1]=a2; atoms[2]=a3; atoms[3]=a4;
}

OutSimpleCo::OutSimpleCo(const Ref<KeyVal> &kv) :
  SimpleCo(kv,4)
{
}

OutSimpleCo::~OutSimpleCo()
{
}

OutSimpleCo&
OutSimpleCo::operator=(const OutSimpleCo& s)
{
  if(label_) delete[] label_;
  label_=new char[strlen(s.label_)+1];
  strcpy(label_,s.label_);
  atoms[0]=s.atoms[0]; atoms[1]=s.atoms[1]; atoms[2]=s.atoms[2];
  atoms[3]=s.atoms[3];
  return *this;
}

double
OutSimpleCo::calc_intco(Molecule& m, double *bmat, double coeff)
{
  int a=atoms[0]-1; int b=atoms[1]-1; int c=atoms[2]-1; int d=atoms[3]-1;
  SCVector3 u1,u2,u3,z1;

  SCVector3 ra(m.r(a));
  SCVector3 rb(m.r(b));
  SCVector3 rc(m.r(c));
  SCVector3 rd(m.r(d));

  u1 = ra-rb;
  u1.normalize();
  u2 = rc-rb;
  u2.normalize();
  u3 = rd-rb;
  u3.normalize();

  z1 = u2.perp_unit(u3);
  double st=u1.dot(z1);
  double ct=s2(st);

  value_ = (st<0) ? -acos(ct) : acos(ct);

  if (bmat) {
    double uu,vv;
    SCVector3 ww,xx,zz;
    double cphi1 = u2.dot(u3);
    double sphi1 = s2(cphi1);
    double cphi2 = u3.dot(u1);
    double cphi3 = u2.dot(u1);
    double den = ct * sphi1*sphi1;
    double sthta2 = (cphi1*cphi2-cphi3)/
              (den*rc.dist(rb));
    double sthta3 = (cphi1*cphi3-cphi2)/
              (den*rd.dist(rb));
#if OLD_BMAT
    sthta2 /= bohr;
    sthta3 /= bohr;
#endif
    int j;
    for(j=0; j < 3; j++) {
      ww[j] = z1[j]*sthta2;
      zz[j] = z1[j]*sthta3;
    }
    xx = z1.perp_unit(u1);
    z1 = u1.perp_unit(xx);
    double r1i = 1.0/ra.dist(rb);
#if OLD_BMAT
    r1i /= bohr;
#endif    
    for(j=0; j < 3; j++) {
      uu = z1[j]*r1i;
      vv = -uu-ww[j]-zz[j];
      bmat[a*3+j] += coeff*uu;
      bmat[b*3+j] += coeff*vv;
      bmat[c*3+j] += coeff*ww[j];
      bmat[d*3+j] += coeff*zz[j];
    }
  }

  return value_;
}


double
OutSimpleCo::calc_force_con(Molecule& m)
{
  int x=atoms[0]-1;
  int a=atoms[1]-1; int b=atoms[2]-1; int c=atoms[3]-1;

  SCVector3 ra(m.r(a));
  SCVector3 rx(m.r(x));

  double rad_ab =   m.atominfo()->atomic_radius(m.Z(a))
                  + m.atominfo()->atomic_radius(m.Z(b));

  double rad_ac =   m.atominfo()->atomic_radius(m.Z(a))
                  + m.atominfo()->atomic_radius(m.Z(c));

  double rad_ax =   m.atominfo()->atomic_radius(m.Z(a))
                  + m.atominfo()->atomic_radius(m.Z(x));

  double r_ax = ra.dist(rx);

  calc_intco(m);

  double k = 0.0025 + 0.0061*pow((rad_ab*rad_ac),0.80)*pow(cos(value()),4.0) *
                           exp(-3.0*(r_ax-rad_ax));

#if OLD_BMAT
  // return force constant in mdyn*ang/rad^2
  return k*4.359813653;
#else  
  return k;
#endif  
}

const char *
OutSimpleCo::ctype() const
{
  return "OUT";
}

double
OutSimpleCo::radians() const
{
  return value_;
}

double
OutSimpleCo::degrees() const
{
  return value_*rtd;
}

double
OutSimpleCo::preferred_value() const
{
  return value_*rtd;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
