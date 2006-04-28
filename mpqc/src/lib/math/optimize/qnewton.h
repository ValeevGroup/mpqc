//
// qnewton.h
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

#ifndef _math_optimize_qnewton_h
#define _math_optimize_qnewton_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/state/state.h>
#include <util/class/class.h>
#include <math/scmat/matrix.h>
#include <math/optimize/function.h>
#include <math/optimize/opt.h>
#include <math/optimize/update.h>

namespace sc {

// //////////////////////////////////////////////////////////////////////
// newton and related methods


/** The QNewtonOpt implements a quasi-Newton optimization scheme. */
class QNewtonOpt: public Optimize {

  protected:
    double maxabs_gradient;
    double accuracy_;

    RefSymmSCMatrix ihessian_;
    Ref<HessianUpdate> update_;
    Ref<LineOpt> lineopt_;

    int take_newton_step_;

    int print_hessian_;
    int print_x_;
    int print_gradient_;
    int linear_;
    int restrict_;
    int dynamic_grad_acc_;
    int force_search_;
    int restart_;

  public:
    /** This KeyVal constructor is used to construct QNewtonOpt
        objects from the input.

        The keywords used by this constructor are listed below.  The KeyVal
        constructor for the parent class, Optimize, will also be
        called, so consult the documentation for
        Optimize(const Ref<KeyVal>&) for additional keywords that
        will be read.

        <dl>

        <dt><tt>update</tt><dd> This gives a HessianUpdate object.  The
        default is to not update the hessian.

        <dt><tt>hessian</tt><dd> By default, the guess hessian is obtained
        from the Function object.  This keyword specifies an lower triangle
        array (the second index must be less than or equal to than the
        first) that replaces the guess hessian.  If some of the elements
        are not given, elements from the guess hessian will be used.

        <dt><tt>lineopt</tt><dd> This gives a LineOpt object for doing line
        optimizations in the Newton direction.  The default is to skip the
        line optimizations.

        <dt><tt>accuracy</tt><dd> The accuracy with which the first
        gradient will be computed.  If this is too large, it may be
        necessary to evaluate the first gradient point twice.  If it is too
        small, it may take longer to evaluate the first point. The default
        is 0.0001.

        <dt><tt>print_x</tt><dd> If true, print the coordinates each
        iteration.  The default is false.

        <dt><tt>print_gradient</tt><dd> If true, print the gradient each
        iteration. The default is false.

        <dt><tt>print_hessian</tt><dd> If true, print the approximate
        hessian each iteration. The default is false.

        <dt><tt>restrict</tt><dd> Use step size restriction when not
        using a line search.  The default is true.

        </dl> */
    QNewtonOpt(const Ref<KeyVal>&);
    QNewtonOpt(StateIn&);
    ~QNewtonOpt();
    void save_data_state(StateOut&);

    void apply_transform(const Ref<NonlinearTransform>&);

    void init();
    int update();
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
