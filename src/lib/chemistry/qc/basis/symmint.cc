//
// symmint.cc
//
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

#include <chemistry/qc/basis/symmint.h>

using namespace sc;

////////////////////////////////////////////////////////////////////////////
// SymmOneBodyIntIter

SymmOneBodyIntIter::SymmOneBodyIntIter(const Ref<OneBodyInt>& ints,
                                       const Ref<PetiteList>& p) :
  OneBodyIntIter(ints), pl(p)
{
}

SymmOneBodyIntIter::~SymmOneBodyIntIter()
{
}

void
SymmOneBodyIntIter::next()
{
  OneBodyIntIter::next();
  while (OneBodyIntIter::ready() && !pl->lambda(icur,jcur))
    OneBodyIntIter::next();
}

void
SymmOneBodyIntIter::start(int ist, int jst, int ien, int jen)
{
  OneBodyIntIter::start(ist,jst,ien,jen);
  while (OneBodyIntIter::ready() && !pl->lambda(icur,jcur))
    OneBodyIntIter::next();
}

double
SymmOneBodyIntIter::scale() const
{
  return (double) pl->lambda(icur,jcur) / (double) pl->order();
}

bool
SymmOneBodyIntIter::cloneable()
{
  return obi->cloneable();
}

Ref<OneBodyIntIter>
SymmOneBodyIntIter::clone()
{
  return new SymmOneBodyIntIter(obi->clone(), pl->clone());
}

////////////////////////////////////////////////////////////////////////////
// SymmTwoBodyIntIter

SymmTwoBodyIntIter::SymmTwoBodyIntIter(const Ref<TwoBodyInt>& ints,
                                       const Ref<PetiteList>& p) :
  TwoBodyIntIter(ints), pl(p)
{
}

SymmTwoBodyIntIter::~SymmTwoBodyIntIter()
{
}

// very inefficient...fix later
void
SymmTwoBodyIntIter::next()
{
  TwoBodyIntIter::next();
  while (TwoBodyIntIter::ready() &&
         !pl->in_p4(i_offset(icur)+jcur, i_offset(kcur)+lcur,
                    icur,jcur,kcur,lcur))
    TwoBodyIntIter::next();
}

// very inefficient...fix later
void
SymmTwoBodyIntIter::start()
{
  TwoBodyIntIter::start();
  while (TwoBodyIntIter::ready() &&
         !pl->in_p4(i_offset(icur)+jcur, i_offset(kcur)+lcur,
                    icur,jcur,kcur,lcur))
    TwoBodyIntIter::next();
}

// very inefficient...fix later
double
SymmTwoBodyIntIter::scale() const
{
  return (double)
    pl->in_p4(i_offset(icur)+jcur, i_offset(kcur)+lcur,icur,jcur,kcur,lcur);
}

////////////////////////////////////////////////////////////////////////////
// SymmTwoBodyTwoCenterIntIter

SymmTwoBodyTwoCenterIntIter::SymmTwoBodyTwoCenterIntIter(const Ref<TwoBodyTwoCenterInt>& ints,
                                       const Ref<PetiteList>& p) :
  TwoBodyTwoCenterIntIter(ints), pl(p)
{
}

SymmTwoBodyTwoCenterIntIter::~SymmTwoBodyTwoCenterIntIter()
{
}

void
SymmTwoBodyTwoCenterIntIter::next()
{
  TwoBodyTwoCenterIntIter::next();
  while (TwoBodyTwoCenterIntIter::ready() && !pl->lambda(icur,jcur))
    TwoBodyTwoCenterIntIter::next();
}

void
SymmTwoBodyTwoCenterIntIter::start(int ist, int jst, int ien, int jen)
{
  TwoBodyTwoCenterIntIter::start(ist,jst,ien,jen);
  while (TwoBodyTwoCenterIntIter::ready() && !pl->lambda(icur,jcur))
    TwoBodyTwoCenterIntIter::next();
}

double
SymmTwoBodyTwoCenterIntIter::scale() const
{
  return (double) pl->lambda(icur,jcur) / (double) pl->order();
}

bool
SymmTwoBodyTwoCenterIntIter::cloneable()
{
  return tbi->cloneable();
}

Ref<TwoBodyTwoCenterIntIter>
SymmTwoBodyTwoCenterIntIter::clone()
{
  return new SymmTwoBodyTwoCenterIntIter(tbi->clone(), pl->clone());
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
