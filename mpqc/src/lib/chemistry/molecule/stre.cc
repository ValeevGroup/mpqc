//
// stre.cc
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

/* stre.cc -- implementation of the stretch internal coordinate class
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

static ClassDesc StreSimpleCo_cd(
  typeid(StreSimpleCo),"StreSimpleCo",1,"public SimpleCo",
  create<StreSimpleCo>, create<StreSimpleCo>, create<StreSimpleCo>);
SimpleCo_IMPL(StreSimpleCo);

StreSimpleCo::StreSimpleCo() : SimpleCo(2) {}

StreSimpleCo::StreSimpleCo(const StreSimpleCo& s)
  : SimpleCo(2)
{
  *this=s;
}

StreSimpleCo::StreSimpleCo(const char *re, int a1, int a2)
  : SimpleCo(2,re)
{
  atoms[0]=a1; atoms[1]=a2;
}

StreSimpleCo::StreSimpleCo(const Ref<KeyVal> &kv)
  : SimpleCo(kv,2)
{
}

StreSimpleCo::~StreSimpleCo()
{
}

StreSimpleCo&
StreSimpleCo::operator=(const StreSimpleCo& s)
{
  if(label_) delete[] label_;
  label_=new char[strlen(s.label_)+1]; strcpy(label_,s.label_);
  atoms[0]=s.atoms[0]; atoms[1]=s.atoms[1];

  return *this;
}

double
StreSimpleCo::calc_force_con(Molecule& m)
{
  int a=atoms[0]-1; int b=atoms[1]-1;
  double rad_ab =   m.atominfo()->atomic_radius(m.Z(a))
                  + m.atominfo()->atomic_radius(m.Z(b));

  calc_intco(m);

  double k = 0.3601 * exp(-1.944*(value()-rad_ab));

#if OLD_BMAT
  // return force constant in mdyn/ang
  return k*4.359813653/(0.52917706*0.52917706);
#else  
  return k;
#endif  
}

double
StreSimpleCo::calc_intco(Molecule& m, double *bmat, double coeff)
{
  int a=atoms[0]-1; int b=atoms[1]-1;
  SCVector3 ra(m.r(a)), rb(m.r(b));
  value_ = ra.dist(rb);
  if(bmat) {
    SCVector3 uu = ra - rb;
    uu.normalize();
    bmat[a*3] += coeff*uu[0]; bmat[b*3] -= coeff*uu[0];
    bmat[a*3+1] += coeff*uu[1]; bmat[b*3+1] -= coeff*uu[1];
    bmat[a*3+2] += coeff*uu[2]; bmat[b*3+2] -= coeff*uu[2];
  }

  return angstrom();
}


const char *
StreSimpleCo::ctype() const
{
  return "STRE";
}

double
StreSimpleCo::bohr() const
{
  return value_;
}

double
StreSimpleCo::angstrom() const
{
  return value_*0.52917706;
}

double
StreSimpleCo::preferred_value() const
{
  return value_*0.52917706;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
