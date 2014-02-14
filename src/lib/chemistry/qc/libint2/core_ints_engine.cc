//
// core_ints_engine.cc
//
// Copyright (C) 2014 Edward Valeev
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

#include<chemistry/qc/libint2/core_ints_engine.h>
#include<chemistry/qc/basis/fjt.h>
#include <libint2/libint2.h>
#include <libint2/boys.h>

using namespace sc;

//typedef sc::FJT _Engine0;
//template<> Ref<CoreIntsEngine<_Engine0>::Engine>
//  CoreIntsEngine<_Engine0>::default_engine_(new CoreIntsEngine<_Engine0>::Engine(12));

typedef ::libint2::FmEval_Chebyshev3 _Engine1;
template<> Ref<CoreIntsEngine<_Engine1>::Engine>
  CoreIntsEngine<_Engine1>::default_engine_(new CoreIntsEngine<_Engine1>::Engine(12));

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
