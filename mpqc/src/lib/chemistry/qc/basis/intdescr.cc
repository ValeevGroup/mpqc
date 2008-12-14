//
// intdescr.cc
//
// Copyright (C) 2005 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
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

#ifdef __GNUG__
#pragma implementation
#endif

#include <util/class/scexception.h>
#include <chemistry/qc/basis/intdescr.h>

using namespace sc;

////

TwoBodyIntDescrERI::TwoBodyIntDescrERI(const Ref<Integral>& IF) :
  TwoBodyIntDescr(), factory_(IF)
{
}

Ref<TwoBodyInt>
TwoBodyIntDescrERI::inteval() const
{
  return factory_->electron_repulsion();
}

unsigned int
TwoBodyIntDescrERI::intset(TwoBodyInt::tbint_type t) const
{
  return TwoBodyIntDescrERI::intSet(t);
}

TwoBodyInt::tbint_type
TwoBodyIntDescrERI::intset(unsigned int t) const
{
  return TwoBodyIntDescrERI::intSet(t);
}

unsigned int
TwoBodyIntDescrERI::intSet(TwoBodyInt::tbint_type t)
{
  switch(t) {
  case TwoBodyInt::eri:
    return 0;
  }
  throw ProgrammingError("TwoBodyIntDescrERI::intSet() -- this type not recognized");
}

TwoBodyInt::tbint_type
TwoBodyIntDescrERI::intSet(unsigned int t)
{
  switch(t) {
  case 0:
    return TwoBodyInt::eri;
  }
  throw ProgrammingError("TwoBodyIntDescrERI::intSet() -- this type not recognized");
}

////

TwoBodyIntDescrR12::TwoBodyIntDescrR12(const Ref<Integral>& IF) :
  TwoBodyIntDescr(), factory_(IF)
{
}

Ref<TwoBodyInt>
TwoBodyIntDescrR12::inteval() const
{
  return factory_->grt();
}

unsigned int
TwoBodyIntDescrR12::intset(TwoBodyInt::tbint_type t) const
{
  return TwoBodyIntDescrR12::intSet(t);
}

TwoBodyInt::tbint_type
TwoBodyIntDescrR12::intset(unsigned int t) const
{
  return TwoBodyIntDescrR12::intSet(t);
}

unsigned int
TwoBodyIntDescrR12::intSet(TwoBodyInt::tbint_type t)
{
    switch(t) {
    case TwoBodyInt::eri:
	return 0;
    case TwoBodyInt::r12:
	return 1;
    case TwoBodyInt::r12t1:
	return 2;
    case TwoBodyInt::r12t2:
	return 3;
    }
  throw ProgrammingError("TwoBodyIntDescrR12::intSet() -- this type not recognized");
}

TwoBodyInt::tbint_type
TwoBodyIntDescrR12::intSet(unsigned int t)
{
    switch(t) {
    case 0:
	return TwoBodyInt::eri;
    case 1:
	return TwoBodyInt::r12;
    case 2:
	return TwoBodyInt::r12t1;
    case 3:
	return TwoBodyInt::r12t2;
    }
  throw ProgrammingError("TwoBodyIntDescrR12::intSet() -- this type not recognized");
}


////

TwoBodyIntDescrG12::TwoBodyIntDescrG12(const Ref<Integral>& IF,
                                         const Ref<IntParamsG12>& params) :
  TwoBodyIntDescr(), factory_(IF), params_(params)
{
}

Ref<TwoBodyInt>
TwoBodyIntDescrG12::inteval() const
{
  return factory_->g12(params_);
}

unsigned int
TwoBodyIntDescrG12::intset(TwoBodyInt::tbint_type t) const
{
  return TwoBodyIntDescrG12::intSet(t);
}

TwoBodyInt::tbint_type
TwoBodyIntDescrG12::intset(unsigned int t) const
{
  return TwoBodyIntDescrG12::intSet(t);
}

unsigned int
TwoBodyIntDescrG12::intSet(TwoBodyInt::tbint_type t)
{
  switch(t) {
  case TwoBodyInt::eri:
    return 0;
  case TwoBodyInt::r12_0_g12:
    return 1;
  case TwoBodyInt::r12_m1_g12:
    return 2;
  case TwoBodyInt::t1g12:
    return 3;
  case TwoBodyInt::t2g12:
    return 4;
  case TwoBodyInt::g12t1g12:
    return 5;
  }
  throw ProgrammingError("TwoBodyIntDescrG12::intSet() -- this type not recognized");
}

TwoBodyInt::tbint_type
TwoBodyIntDescrG12::intSet(unsigned int t)
{
  switch(t) {
  case 0:
    return TwoBodyInt::eri;
  case 1:
    return TwoBodyInt::r12_0_g12;
  case 2:
    return TwoBodyInt::r12_m1_g12;
  case 3:
    return TwoBodyInt::t1g12;
  case 4:
    return TwoBodyInt::t2g12;
  case 5:
    return TwoBodyInt::g12t1g12;
  }
  throw ProgrammingError("TwoBodyIntDescrG12::intSet() -- this type not recognized");
}

////

TwoBodyIntDescrG12NC::TwoBodyIntDescrG12NC(const Ref<Integral>& IF,
					   const Ref<IntParamsG12>& params) :
  TwoBodyIntDescr(), factory_(IF), params_(params)
{
}

Ref<TwoBodyInt>
TwoBodyIntDescrG12NC::inteval() const
{
  return factory_->g12nc(params_);
}

unsigned int
TwoBodyIntDescrG12NC::intset(TwoBodyInt::tbint_type t) const
{
  return TwoBodyIntDescrG12NC::intSet(t);
}

TwoBodyInt::tbint_type
TwoBodyIntDescrG12NC::intset(unsigned int t) const
{
  return TwoBodyIntDescrG12NC::intSet(t);
}

unsigned int
TwoBodyIntDescrG12NC::intSet(TwoBodyInt::tbint_type t)
{
  switch(t) {
  case TwoBodyInt::eri:
    return 0;
  case TwoBodyInt::r12_0_g12:
    return 1;
  case TwoBodyInt::r12_m1_g12:
    return 2;
  case TwoBodyInt::g12t1g12:
    return 3;
  case TwoBodyInt::anti_g12g12:
    return 4;
  }
  throw ProgrammingError("TwoBodyIntDescrG12NC::intSet() -- this type not recognized");
}

TwoBodyInt::tbint_type
TwoBodyIntDescrG12NC::intSet(unsigned int t)
{
  switch(t) {
  case 0:
    return TwoBodyInt::eri;
  case 1:
    return TwoBodyInt::r12_0_g12;
  case 2:
    return TwoBodyInt::r12_m1_g12;
  case 3:
    return TwoBodyInt::g12t1g12;
  case 4:
    return TwoBodyInt::anti_g12g12;
  }
  throw ProgrammingError("TwoBodyIntDescrG12NC::intSet() -- this type not recognized");
}

////

TwoBodyIntDescrG12DKH::TwoBodyIntDescrG12DKH(const Ref<Integral>& IF,
                                             const Ref<IntParamsG12>& params) :
  TwoBodyIntDescr(), factory_(IF), params_(params)
{
  // currently implemented only when bra and ket functions are the same
  if (params_->bra() != params_->ket())
    throw FeatureNotImplemented("TwoBodyIntDescrG12DKH -- g12dkh integrals only implemented when only 1 correlation factor is used");
}

Ref<TwoBodyInt>
TwoBodyIntDescrG12DKH::inteval() const
{
  return factory_->g12dkh(params_);
}

unsigned int
TwoBodyIntDescrG12DKH::intset(TwoBodyInt::tbint_type t) const
{
  return TwoBodyIntDescrG12DKH::intSet(t);
}

TwoBodyInt::tbint_type
TwoBodyIntDescrG12DKH::intset(unsigned int t) const
{
  return TwoBodyIntDescrG12DKH::intSet(t);
}

unsigned int
TwoBodyIntDescrG12DKH::intSet(TwoBodyInt::tbint_type t)
{
  switch(t) {
  case TwoBodyInt::g12p4g12_m_g12t1g12t1:
    return 0;
  }
  throw ProgrammingError("TwoBodyIntDescrG12DKH::intSet() -- this type not recognized");
}

TwoBodyInt::tbint_type
TwoBodyIntDescrG12DKH::intSet(unsigned int t)
{
  switch(t) {
  case 0:
    return TwoBodyInt::g12p4g12_m_g12t1g12t1;
  }
  throw ProgrammingError("TwoBodyIntDescrG12DKH::intSet() -- this type not recognized");
}

////

TwoBodyIntDescrGenG12::TwoBodyIntDescrGenG12(const Ref<Integral>& IF,
					     const Ref<IntParamsGenG12>& params) :
  TwoBodyIntDescr(), factory_(IF), params_(params)
{
}

Ref<TwoBodyInt>
TwoBodyIntDescrGenG12::inteval() const
{
  return factory_->geng12(params_);
}

unsigned int
TwoBodyIntDescrGenG12::intset(TwoBodyInt::tbint_type t) const
{
  return TwoBodyIntDescrGenG12::intSet(t);
}

TwoBodyInt::tbint_type
TwoBodyIntDescrGenG12::intset(unsigned int t) const
{
  return TwoBodyIntDescrGenG12::intSet(t);
}

unsigned int
TwoBodyIntDescrGenG12::intSet(TwoBodyInt::tbint_type t)
{
  switch(t) {
  case TwoBodyInt::eri:
    return 0;
  case TwoBodyInt::r12_0_gg12:
    return 1;
  case TwoBodyInt::r12_m1_gg12:
    return 2;
  case TwoBodyInt::gg12t1gg12:
    return 3;
  }
  throw ProgrammingError("TwoBodyIntDescrGenG12::intSet() -- this type not recognized");
}

TwoBodyInt::tbint_type
TwoBodyIntDescrGenG12::intSet(unsigned int t)
{
  switch(t) {
  case 0:
    return TwoBodyInt::eri;
  case 1:
    return TwoBodyInt::r12_0_gg12;
  case 2:
    return TwoBodyInt::r12_m1_gg12;
  case 3:
    return TwoBodyInt::gg12t1gg12;
  }
  throw ProgrammingError("TwoBodyIntDescrGenG12::intSet() -- this type not recognized");
}

