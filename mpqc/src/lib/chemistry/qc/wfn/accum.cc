//
// accum.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <chemistry/qc/wfn/accum.h>
#include <chemistry/qc/basis/integral.h>

///////////////////////////////////////////////////////////////////////////
// AccumDIH

#define CLASSNAME AccumDIH
#define PARENTS public SavableState
#include <util/class/classia.h>

void *
AccumDIH::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

AccumDIH::AccumDIH()
{
}

AccumDIH::AccumDIH(StateIn&s) :
  SavableState(s)
{
  basis_set_.restore_state(s);
  integral_.restore_state(s);
}

// for now I'm assuming that the specializations of AccumDIH will call
// init() so we won't read in the integral or basis set here
AccumDIH::AccumDIH(const RefKeyVal&)
{
}

AccumDIH::~AccumDIH()
{
}

void
AccumDIH::save_data_state(StateOut& s)
{
  basis_set_.save_state(s);
  integral_.save_state(s);
}

void
AccumDIH::init(const RefGaussianBasisSet& b, const RefIntegral& i)
{
  basis_set_ = b;
  integral_ = i;
}

void
AccumDIH::done()
{
}

///////////////////////////////////////////////////////////////////////////
// AccumDDH

#define CLASSNAME AccumDDH
#define PARENTS public SavableState
#include <util/class/classia.h>
void *
AccumDDH::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

AccumDDH::AccumDDH()
{
}

AccumDDH::AccumDDH(StateIn& s) :
  SavableState(s)
{
  basis_set_.restore_state(s);
  integral_.restore_state(s);
}

// for now I'm assuming that the specializations of AccumDDH will call
// init() so we won't read in the integral or basis set here
AccumDDH::AccumDDH(const RefKeyVal&)
{
}

AccumDDH::~AccumDDH()
{
}

void
AccumDDH::save_data_state(StateOut& s)
{
  basis_set_.save_state(s);
  integral_.save_state(s);
}

void
AccumDDH::init(const RefGaussianBasisSet& b, const RefIntegral& i)
{
  basis_set_ = b;
  integral_ = i;
}

void
AccumDDH::done()
{
}

///////////////////////////////////////////////////////////////////////////
// AccumDDH

#define CLASSNAME AccumNullDDH
#define PARENTS public AccumDDH
#define HAVE_CTOR
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>

void *
AccumNullDDH::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = AccumDDH::_castdown(cd);
  return do_castdowns(casts,cd);
}

AccumNullDDH::AccumNullDDH()
{
}

AccumNullDDH::AccumNullDDH(StateIn&s) :
  AccumDDH(s)
{
}

AccumNullDDH::AccumNullDDH(const RefKeyVal& keyval) :
  AccumDDH(keyval)
{
}

AccumNullDDH::~AccumNullDDH()
{
}

void
AccumNullDDH::save_data_state(StateOut& s)
{
  AccumDDH::save_data_state(s);
}

void
AccumNullDDH::accum(const RefSymmSCMatrix& h, const RefSymmSCMatrix& h_open)
{
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")
// End:
