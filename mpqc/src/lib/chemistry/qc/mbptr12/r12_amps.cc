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

const RefSCMatrix R12Amplitudes::T2_aa() const  {return T2_aa_;};
const RefSCMatrix R12Amplitudes::T2_ab() const  {return T2_ab_;};
const RefSCMatrix R12Amplitudes::Rvv_aa() const {return Rvv_aa_;};
const RefSCMatrix R12Amplitudes::Rvv_ab() const {return Rvv_ab_;};
const RefSCMatrix R12Amplitudes::Roo_aa() const {return Roo_aa_;};
const RefSCMatrix R12Amplitudes::Roo_ab() const {return Roo_ab_;};
const RefSCMatrix R12Amplitudes::Rvo_aa() const {return Rvo_aa_;};
const RefSCMatrix R12Amplitudes::Rvo_ab() const {return Rvo_ab_;};
const RefSCMatrix R12Amplitudes::Rxo_aa() const {return Rxo_aa_;};
const RefSCMatrix R12Amplitudes::Rxo_ab() const {return Rxo_ab_;};




/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
