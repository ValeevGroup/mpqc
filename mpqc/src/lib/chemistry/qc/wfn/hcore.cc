//
// hcore.cc
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <chemistry/qc/wfn/hcore.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/symmint.h>

#define CLASSNAME AccumHCore
#define PARENTS public AccumDIH
#define HAVE_CTOR
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>

void *
AccumHCore::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = AccumDIH::_castdown(cd);
  return do_castdowns(casts,cd);
}

AccumHCore::AccumHCore()
{
}

AccumHCore::AccumHCore(StateIn&s) :
  AccumDIH(s)
{
}

AccumHCore::AccumHCore(const RefKeyVal& keyval) :
  AccumDIH(keyval)
{
}

AccumHCore::~AccumHCore()
{
}

void
AccumHCore::save_data_state(StateOut& s)
{
  AccumDIH::save_data_state(s);
}

void
AccumHCore::accum(const RefSymmSCMatrix& h)
{
  integral_->set_basis(basis_set_);

  RefSCElementOp hc = new OneBodyIntOp(integral_->kinetic());
  h.assign(0.0);
  h.element_op(hc);
  hc=0;

  RefOneBodyInt nuc = integral_->nuclear();
  nuc->reinitialize();
  hc = new OneBodyIntOp(nuc);
  h.element_op(hc);
  hc=0;
}

//////////////////////////////////////////////////////////////////////////////

#define CLASSNAME SymmAccumHCore
#define PARENTS public AccumDIH
#define HAVE_CTOR
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>

void *
SymmAccumHCore::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = AccumDIH::_castdown(cd);
  return do_castdowns(casts,cd);
}

SymmAccumHCore::SymmAccumHCore()
{
}

SymmAccumHCore::SymmAccumHCore(StateIn&s) :
  AccumDIH(s)
{
}

SymmAccumHCore::SymmAccumHCore(const RefKeyVal& keyval) :
  AccumDIH(keyval)
{
}

SymmAccumHCore::~SymmAccumHCore()
{
}

void
SymmAccumHCore::save_data_state(StateOut& s)
{
  AccumDIH::save_data_state(s);
}

void
SymmAccumHCore::accum(const RefSymmSCMatrix& h)
{
  integral_->set_basis(basis_set_);
  RefPetiteList pl = integral_->petite_list();

  // form skeleton Hcore in AO basis
  RefSymmSCMatrix hao(basis_set_->basisdim(), basis_set_->matrixkit());
  hao.assign(0.0);

  RefSCElementOp hc =
    new OneBodyIntOp(new SymmOneBodyIntIter(integral_->kinetic(), pl));
  hao.element_op(hc);
  hc=0;

  RefOneBodyInt nuc = integral_->nuclear();
  nuc->reinitialize();
  hc = new OneBodyIntOp(new SymmOneBodyIntIter(nuc, pl));
  hao.element_op(hc);
  hc=0;

  // now symmetrize Hao
  pl->symmetrize(hao,h);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
