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
#include <math/optimize/transform.h>

namespace sc {

// //////////////////////////////////////////////////////////////////////
//  hessian update classes.  based on the value of inverse_hessian_
//  x and g may be reversed (see Schlegel, ab initio Methods in Quantum
//  Chemistry I, 1987, p 10


/** The HessianUpdate abstract class is used to specify a hessian update
scheme.  It is used, for example, by QNewtonOpt objects.  Based on the
value of inverse_hessian_ x and g may be reversed (see Schlegel, Ab initio
Methods in Quantum Chemistry I, 1987, p 10).
*/
class HessianUpdate: virtual public SavableState {
  protected:
    int inverse_hessian_;
  public:
    HessianUpdate();
    HessianUpdate(StateIn&);
    HessianUpdate(const Ref<KeyVal>&);
    void save_data_state(StateOut&);
    virtual ~HessianUpdate();
    virtual void update(const RefSymmSCMatrix&hessian,const Ref<Function>&,
                        const RefSCVector&xnew,const RefSCVector&gnew) = 0;
    virtual void set_inverse();
    virtual void apply_transform(const Ref<NonlinearTransform>&);
};


/** The DFPUpdate class is used to specify a Davidson, Fletcher, and Powell
hessian update scheme.  */
class DFPUpdate: public HessianUpdate {
  protected:
    RefSCVector xprev;
    RefSCVector gprev;
  public:
    DFPUpdate();
    DFPUpdate(StateIn&);
    /** The KeyVal constructor reads the following keywords:
    
    <dl>

    <dt><tt>xprev</tt><dd> The previous coordinates can be given (but is not
    recommended).  The default is none.

    <dt><tt>gprev</tt><dd> The previous gradient can be given (but is not
    recommended).  The default is none.

    </dl>

    */
    DFPUpdate(const Ref<KeyVal>&);
    void save_data_state(StateOut&);
    ~DFPUpdate();
    void update(const RefSymmSCMatrix&ihessian,const Ref<Function>&,
                const RefSCVector&xnew,const RefSCVector&gnew);
    void apply_transform(const Ref<NonlinearTransform>&);
    void set_inverse();
};

/** The DFPUpdate class is used to specify a Broyden, Fletcher, Goldfarb,
and Shanno hessian update scheme.  This hessian update method is the
recommended method for use with QNewtonOpt objects.  */
class BFGSUpdate: public DFPUpdate {
  public:
    BFGSUpdate();
    BFGSUpdate(StateIn&);
    BFGSUpdate(const Ref<KeyVal>&);
    void save_data_state(StateOut&);
    ~BFGSUpdate();
    void update(const RefSymmSCMatrix&ihessian,const Ref<Function>&,
                const RefSCVector&xnew,const RefSCVector&gnew);
};

/** The PowellUpdate class is used to specify a Powell hessian update.
This hessian update method is the recommended method for use with
transition state searches (the EFCOpt class can be used for transition
state searches).  */
class PowellUpdate: public HessianUpdate {
  protected:
    RefSCVector xprev;
    RefSCVector gprev;
  public:
    PowellUpdate();
    PowellUpdate(StateIn&);
    PowellUpdate(const Ref<KeyVal>&);
    void save_data_state(StateOut&);
    ~PowellUpdate();
    void update(const RefSymmSCMatrix&ihessian,const Ref<Function>&func,
                const RefSCVector&xnew,const RefSCVector&gnew);
    void apply_transform(const Ref<NonlinearTransform>&);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
