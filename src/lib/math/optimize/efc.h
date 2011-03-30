//
// efc.h
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

#ifndef _math_optimize_efc_h
#define _math_optimize_efc_h

#include <util/state/state.h>
#include <util/class/class.h>
#include <math/scmat/matrix.h>
#include <math/optimize/function.h>
#include <math/optimize/opt.h>
#include <math/optimize/update.h>

namespace sc {

/**
The EFCOpt class implements
eigenvector following as described by Baker in J. Comput. Chem., Vol 7, No
4, 385-395, 1986.
*/
class EFCOpt: public Optimize {
  protected:
    int tstate;
    int modef;
  
    double maxabs_gradient;
    double convergence_;
    double accuracy_;

    RefSymmSCMatrix hessian_;
    Ref<HessianUpdate> update_;
    RefSCVector last_mode_;

  public:
    /** The KeyVal constructor reads the following keywords:

        <dl>

        <dt><tt>update</tt><dd> This gives an HessianUpdate object.  The
        default is to not update the hessian.

        <dt><tt>transition_state</tt><dd> If this is true than a transition
        state search will be performed. The default is false.

        <dt><tt>mode_following</tt><dd> If this is true, then the initial
        search direction for a transition state search will be choosen to
        similar to the first coordinate of the Function.  The default is
        false.

        <dt><tt>hessian</tt><dd> By default, the guess hessian is obtained
        from the Function object.  This keyword specifies an lower triangle
        array (the second index must be less than or equal to than the
        first) that replaces the guess hessian.  If some of the elements
        are not given, elements from the guess hessian will be used.

        <dt><tt>accuracy</tt><dd> The accuracy with which the first
        gradient will be computed.  If this is too large, it may be
        necessary to evaluate the first gradient point twice.  If it is too
        small, it may take longer to evaluate the first point. The default
        is 0.0001.

        </dl> */
    EFCOpt(const Ref<KeyVal>&);
    EFCOpt(StateIn&);
    ~EFCOpt();
    void save_data_state(StateOut&);

    void apply_transform(const Ref<NonlinearTransform>&);

    void init();
    int update();

    void print(std::ostream& = ExEnv::out0()) const;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
