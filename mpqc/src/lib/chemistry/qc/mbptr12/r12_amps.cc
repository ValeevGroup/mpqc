//
// r12_amps.cc
//
// Copyright (C) 2004 Edward Valeev
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

#include <chemistry/qc/mbptr12/r12_amps.h>

using namespace std;
using namespace sc;

RefSCMatrix
F12Amplitudes::T2(SpinCase2 S) const
{
  return T2_[S];
}

RefSCMatrix
F12Amplitudes::Fvv(SpinCase2 S) const
{
  return Fvv_[S];
}

RefSCMatrix
F12Amplitudes::Foo(SpinCase2 S) const
{
  return Foo_[S];
}

RefSCMatrix
F12Amplitudes::Fvo(SpinCase2 S) const
{
  return Fvo_[S];
}

RefSCMatrix
F12Amplitudes::Fxo(SpinCase2 S) const
{
  return Fxo_[S];
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
