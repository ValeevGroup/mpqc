//
// accum.h
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

#ifndef _chemistry_qc_wfn_accum_h
#define _chemistry_qc_wfn_accum_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/wfn/wfn.h>

// //////////////////////////////////////////////////////////////////////////

// computes additions to H
class AccumH: virtual public SavableState {
#   define CLASSNAME AccumH
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    RefWavefunction wfn_;

  public:
    AccumH();
    AccumH(StateIn&);
    AccumH(const RefKeyVal&);
    virtual ~AccumH();

    void save_data_state(StateOut&);
    
    virtual void init(const RefWavefunction&);
    virtual void accum(const RefSymmSCMatrix& h) =0;
    virtual void print_summary();
    virtual void done();

    // Returns the scalar contribution to the energy.
    // Available only after accum is called.
    virtual double e();
};
SavableState_REF_dec(AccumH);

class AccumHNull: public AccumH {
#   define CLASSNAME AccumHNull
#   define HAVE_CTOR
#   define HAVE_STATEIN_CTOR
#   define HAVE_KEYVAL_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    AccumHNull();
    AccumHNull(StateIn&);
    AccumHNull(const RefKeyVal&);
    ~AccumHNull();

    void save_data_state(StateOut&);
    
    void accum(const RefSymmSCMatrix& h);
};

class SumAccumH: public AccumH {
#   define CLASSNAME SumAccumH
#   define HAVE_STATEIN_CTOR
#   define HAVE_KEYVAL_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    int n_;
    RefAccumH *accums_;

  public:
    SumAccumH(StateIn&);
    SumAccumH(const RefKeyVal&);
    ~SumAccumH();

    void save_data_state(StateOut&);

    void init(const RefWavefunction&);
    void accum(const RefSymmSCMatrix& h);
    void done();

    double e();
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
