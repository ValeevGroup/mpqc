//
// update.h
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

#ifndef _math_optimize_update_h
#define _math_optimize_update_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/state/state.h>
#include <util/class/class.h>
#include <math/scmat/matrix.h>
#include <math/optimize/function.h>

class RefNonlinearTransform;

////////////////////////////////////////////////////////////////////////
//  hessian update classes.  based on the value of inverse_hessian_
//  x and g may be reversed (see Schlegel, ab initio Methods in Quantum
//  Chemistry I, 1987, p 10

class HessianUpdate: virtual_base public SavableState {
#   define CLASSNAME HessianUpdate
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    int inverse_hessian_;
  public:
    HessianUpdate();
    HessianUpdate(StateIn&);
    HessianUpdate(const RefKeyVal&);
    void save_data_state(StateOut&);
    virtual ~HessianUpdate();
    virtual void update(const RefSymmSCMatrix&hessian,const RefFunction&,
                        const RefSCVector&xnew,const RefSCVector&gnew) = 0;
    virtual void set_inverse();
    virtual void apply_transform(const RefNonlinearTransform&);
};
SavableState_REF_dec(HessianUpdate);

class DFPUpdate: public HessianUpdate {
#   define CLASSNAME DFPUpdate
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    RefSCVector xprev;
    RefSCVector gprev;
  public:
    DFPUpdate();
    DFPUpdate(StateIn&);
    DFPUpdate(const RefKeyVal&);
    void save_data_state(StateOut&);
    ~DFPUpdate();
    void update(const RefSymmSCMatrix&ihessian,const RefFunction&,
                const RefSCVector&xnew,const RefSCVector&gnew);
    void apply_transform(const RefNonlinearTransform&);
    void set_inverse();
};

class BFGSUpdate: public DFPUpdate {
#   define CLASSNAME BFGSUpdate
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    BFGSUpdate();
    BFGSUpdate(StateIn&);
    BFGSUpdate(const RefKeyVal&);
    void save_data_state(StateOut&);
    ~BFGSUpdate();
    void update(const RefSymmSCMatrix&ihessian,const RefFunction&,
                const RefSCVector&xnew,const RefSCVector&gnew);
};

class PowellUpdate: public HessianUpdate {
#   define CLASSNAME PowellUpdate
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    RefSCVector xprev;
    RefSCVector gprev;
  public:
    PowellUpdate();
    PowellUpdate(StateIn&);
    PowellUpdate(const RefKeyVal&);
    void save_data_state(StateOut&);
    ~PowellUpdate();
    void update(const RefSymmSCMatrix&ihessian,const RefFunction&func,
                const RefSCVector&xnew,const RefSCVector&gnew);
    void apply_transform(const RefNonlinearTransform&);
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
