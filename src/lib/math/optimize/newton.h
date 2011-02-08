//
// newton.h
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

#ifndef _math_optimize_newton_h
#define _math_optimize_newton_h

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

/// Implements Newton method
class NewtonOpt: public Optimize {
  protected:
    double maxabs_gradient;
    double accuracy_;

    int print_hessian_;
    int print_x_;
    int print_gradient_;
  public:
    ///
    NewtonOpt(const Ref<KeyVal>&);
    NewtonOpt(StateIn&);
    ~NewtonOpt();
    void save_data_state(StateOut&);

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
