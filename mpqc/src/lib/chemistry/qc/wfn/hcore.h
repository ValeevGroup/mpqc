//
// hcore.h
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

#ifndef _chemistry_qc_wfn_hcore_h
#define _chemistry_qc_wfn_hcore_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/wfn/accum.h>

class AccumHCore: public AccumDIH {
#   define CLASSNAME AccumHCore
#   define HAVE_CTOR
#   define HAVE_STATEIN_CTOR
#   define HAVE_KEYVAL_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    AccumHCore();
    AccumHCore(StateIn&);
    AccumHCore(const RefKeyVal&);
    ~AccumHCore();

    void save_data_state(StateOut&);
    
    void accum(const RefSymmSCMatrix& h);
};
SavableState_REF_dec(AccumHCore);

class SymmAccumHCore: public AccumDIH {
#   define CLASSNAME SymmAccumHCore
#   define HAVE_CTOR
#   define HAVE_STATEIN_CTOR
#   define HAVE_KEYVAL_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    SymmAccumHCore();
    SymmAccumHCore(StateIn&);
    SymmAccumHCore(const RefKeyVal&);
    ~SymmAccumHCore();

    void save_data_state(StateOut&);
    
    void accum(const RefSymmSCMatrix& h);
};
SavableState_REF_dec(SymmAccumHCore);

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")
// End:
