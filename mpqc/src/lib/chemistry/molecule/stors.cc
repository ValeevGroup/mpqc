//
// stors.cc
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

#include <string.h>
#include <util/misc/math.h>

#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <chemistry/molecule/simple.h>
#include <chemistry/molecule/localdef.h>

using namespace sc;

static ClassDesc ScaledTorsSimpleCo_cd(
  typeid(ScaledTorsSimpleCo),"ScaledTorsSimpleCo",1,"public SimpleCo",
  create<ScaledTorsSimpleCo>, create<ScaledTorsSimpleCo>, create<ScaledTorsSimpleCo>);
SimpleCo_IMPL_eq(ScaledTorsSimpleCo)

ScaledTorsSimpleCo::ScaledTorsSimpleCo(StateIn&s):
  SimpleCo(s)
{
  s.get(old_torsion_);
}

void
ScaledTorsSimpleCo::save_data_state(StateOut&s)
{
  SimpleCo::save_data_state(s);
  s.put(old_torsion_);
}

ScaledTorsSimpleCo::ScaledTorsSimpleCo() : SimpleCo(4)
{
  old_torsion_ = 0.0;
}

ScaledTorsSimpleCo::ScaledTorsSimpleCo(const ScaledTorsSimpleCo& s)
  : SimpleCo(4)
{
  *this=s;
  old_torsion_ = 0.0;
}

ScaledTorsSimpleCo::ScaledTorsSimpleCo(const char *refr,
                                       int a1, int a2, int a3, int a4)
  : SimpleCo(4,refr)
{
  atoms[0]=a1; atoms[1]=a2; atoms[2]=a3; atoms[3]=a4;
  old_torsion_ = 0.0;
}

ScaledTorsSimpleCo::~ScaledTorsSimpleCo()
{
}

ScaledTorsSimpleCo::ScaledTorsSimpleCo(const Ref<KeyVal> &kv):
  SimpleCo(kv,4)
{
  old_torsion_ = 0.0;
}

ScaledTorsSimpleCo&
ScaledTorsSimpleCo::operator=(const ScaledTorsSimpleCo& s)
{
  if(label_) delete[] label_;
  label_=new char[strlen(s.label_)+1];
  strcpy(label_,s.label_);
  atoms[0]=s.atoms[0]; atoms[1]=s.atoms[1]; atoms[2]=s.atoms[2];
  atoms[3]=s.atoms[3];
  return *this;
}

double
ScaledTorsSimpleCo::calc_intco(Molecule& m, double *bmat, double coeff)
{
  int a=atoms[0]-1; int b=atoms[1]-1; int c=atoms[2]-1; int d=atoms[3]-1;
  SCVector3 u1,u2,u3,z1,z2;

  SCVector3 ra(m.r(a)), rb(m.r(b)), rc(m.r(c)), rd(m.r(d));

  double rab, rbc, rcd;
  u1 = ra - rb; rab = u1.norm(); u1 *= 1.0/rab;
  u2 = rc - rb; rbc = u2.norm(); u2 *= 1.0/rbc;
  u3 = rc - rd; rcd = u3.norm(); u3 *= 1.0/rcd;

  z1 = u1.perp_unit(u2);
  z2 = u3.perp_unit(u2);

  double co=z1.dot(z2);
  SCVector3 z1xz2 = z1.cross(z2);
  double co2=z1xz2.dot(u2);

  if (co < -1.0) co= -1.0;
  if (co > 1.0) co = 1.0;

  double tors_value = (co2<0) ? acos(-co) : -acos(-co);

  // ok, we want omega between 3*pi/2 and -pi/2, so if omega is > pi/2
  // (omega is eventually -omega), then knock 2pi off of it
  if(tors_value > pih) tors_value -= tpi;

  // the following tests to see if the new coordinate has crossed the
  // 3pi/2 <--> -pi/2 boundary...if so, then we add or subtract 2pi as
  // needed to prevent the transformation from internals to cartesians
  // from blowing up
  while(old_torsion_ - tors_value > pi + 1.0e-6) tors_value += tpi;
  while(old_torsion_ - tors_value < -(pi + 1.0e-6)) tors_value -= tpi;

  // This differs from a normal torsion by the factor
  // scalar(u1,u2)*scalar(u2,u3).  This prevents the
  // value from being ill defined at nearly linear geometries.
  double cos_abc = u1.dot(u2);
  double cos_bcd = u2.dot(u3);
  double sin_abc=s2(cos_abc);
  double sin_bcd=s2(cos_bcd);
  double colin = sin_abc * sin_bcd;
  value_ = tors_value * colin;

  if (bmat) {
    double uu,vv,ww,zz;
    double r1 = rab;
    double r2 = rbc;
    double r3 = rcd;
#if OLD_BMAT
    r1 *= bohr;
    r2 *= bohr;
    r3 *= bohr;
#endif    
    for (int j=0; j < 3; j++) {
      // compute the derivatives for a normal torsion
      if (sin_abc > 1.0e-5) uu = z1[j]/(r1*sin_abc);
      else  uu = 0.0;
      if (sin_bcd > 1.0e-5) zz = z2[j]/(r3*sin_bcd);
      else zz = 0.0;
      vv = (r1*cos_abc/r2-1.0)*uu-zz*r3*cos_bcd/r2;
      // use the chain rule to get the values for the scaled torsion
      static int first = 0;
      if (first) {
          ExEnv::out0() << indent
               << scprintf("uu = %12.8f colin = %12.8f sin_abc = %12.8f\n",
                           uu, colin, sin_abc)
               << indent
               << scprintf("tors_value = %12.8f (%12.8f deg)\n", tors_value,
                           tors_value * 180.0/M_PI)
               << indent
               << scprintf("cos_abc = %12.8f cos_bcd = %12.8f\n",
                           cos_abc, cos_bcd);
        }
      uu = uu*colin;
      if (sin_abc > 1.0e-5) uu += tors_value
                                * (-cos_abc/sin_abc)
                                * sin_bcd
                                * (u2[j] - cos_abc * u1[j])/rab;
      vv = vv*colin;
      if (sin_abc > 1.0e-5) vv += tors_value
                                * (-cos_abc/sin_abc)
                                * sin_bcd
                                * ((-u2[j] + cos_abc*u1[j])/rab
                                   +(-u1[j] + cos_abc*u2[j])/rbc);
      if (sin_bcd > 1.0e-5) vv += tors_value
                                * (-cos_bcd/sin_bcd)
                                * sin_abc
                                * (-u3[j] + cos_bcd * u2[j])/rbc;
      zz = zz*colin;
      if (sin_bcd > 1.0e-5) zz += tors_value
                                * (-cos_bcd/sin_bcd)
                                * sin_abc
                                * (-u2[j] + cos_bcd * u3[j])/rcd;
      ww = -uu-vv-zz;
      bmat[a*3+j] += coeff*uu;
      bmat[b*3+j] += coeff*vv;
      bmat[c*3+j] += coeff*ww;
      bmat[d*3+j] += coeff*zz;
      if (first) first = 0;
    }
  }

  // save the old value of the torsion so we can make sure the discontinuity
  // at -pi/2 doesn't bite us
  old_torsion_ = tors_value;
  return value_;
}

double
ScaledTorsSimpleCo::calc_force_con(Molecule& m)
{
  int a=atoms[1]-1; int b=atoms[2]-1;

  SCVector3 ra(m.r(a));
  SCVector3 rb(m.r(b));

  double rad_ab =   m.atominfo()->atomic_radius(m.Z(a))
                  + m.atominfo()->atomic_radius(m.Z(b));

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
ScaledTorsSimpleCo::ctype() const
{
  return "STOR";
}

double
ScaledTorsSimpleCo::radians() const
{
  return value_;
}

double
ScaledTorsSimpleCo::degrees() const
{
  return value_*rtd;
}

double
ScaledTorsSimpleCo::preferred_value() const
{
  return value_*rtd;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
