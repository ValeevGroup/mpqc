//
// conv.h
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

#ifndef _math_optimize_conv_h
#define _math_optimize_conv_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/state/state.h>
#include <util/class/class.h>
#include <math/scmat/matrix.h>
#include <math/optimize/function.h>

// //////////////////////////////////////////////////////////////////////

/** The Convergence class is used to check for the
    convergence of an optimization. */
class Convergence: virtual_base public SavableState {
#   define CLASSNAME Convergence
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    RefSCVector grad_;
    RefSCVector x_;
    RefSCVector nextx_;
    int use_max_disp_;
    double max_disp_;
    int use_max_grad_;
    double max_grad_;
    int use_rms_disp_;
    double rms_disp_;
    int use_rms_grad_;
    double rms_grad_;
    int use_graddisp_;
    double graddisp_;

    void check_conv(const char *heading,
                    double val, double bound,
                    int &pass, int &fail);

    void set_defaults();
  public:
    /// Standard constructors and destructor.
    Convergence();
    Convergence(double tolerance);
    Convergence(StateIn&);
    Convergence(const RefKeyVal&);
    virtual ~Convergence();

    void save_data_state(StateOut&);

    /// Set the current gradient and displacement.
    virtual void get_grad(const RefFunction &);
    virtual void get_x(const RefFunction &);
    virtual void set_nextx(const RefSCVector &);

    /// Set the current gradient and displacement to null.
    virtual void reset();

    /// Return nonzero if the optimization has converged.
    virtual int converged();
};
SavableState_REF_dec(Convergence);

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
