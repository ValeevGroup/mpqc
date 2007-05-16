//
// intdescr.cc
//
// Copyright (C) 2005 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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


