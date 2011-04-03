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

#include <util/state/state.h>
#include <util/class/class.h>
#include <math/scmat/matrix.h>
#include <math/optimize/function.h>

namespace sc {

// //////////////////////////////////////////////////////////////////////

/** The Convergence class is used by the optimizer to determine when an
optimization is converged.  The KeyVal input for Convergence is given
below.  Giving none of these keywords is the same as giving the following
input:
<pre>
  conv<Convergence>: (
    max_disp = 1.0e-6
    max_grad = 1.0e-6
    graddisp = 1.0e-6
  )
</pre>
*/
class Convergence: virtual public SavableState {
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
    Convergence();
    Convergence(double tolerance);
    Convergence(StateIn&);

    /** The KeyVal constructor reads the following keywords:

        <dl>

        <dt><tt>max_disp</tt><dd> The value of the maximum displacement
        must be less then the value of this keyword for the calculation to
        be converged.  The default is to not check this parameter.
        However, if no other keyword are given, default convergence
        parameters are chosen as described above.

        <dt><tt>max_grad</tt><dd> The value of the maximum gradient must be
        less then the value of this keyword for the calculation to be
        converged.  The default is to not check this parameter.  However,
        if no other keyword are given, default convergence parameters are
        chosen as described above.

        <dt><tt>rms_disp</tt><dd> The value of the RMS of the displacements
        must be less then the value of this keyword for the calculation to
        be converged.  The default is to not check this parameter.
        However, if no other keyword are given, default convergence
        parameters are chosen as described above.

        <dt><tt>rms_grad</tt><dd> The value of the RMS of the gradients
        must be less then the value of this keyword for the calculation to
        be converged.  The default is to not check this parameter.
        However, if no other keyword are given, default convergence
        parameters are chosen as described above.

        <dt><tt>graddisp</tt><dd> The value of the scalar product of the
        gradient vector with the displacement vector must be less then the
        value of this keyword for the calculation to be converged.  The
        default is to not check this parameter.  However, if no other
        keyword are given, default convergence parameters are chosen as
        described above.

        </dl> */
    Convergence(const Ref<KeyVal>&);
    virtual ~Convergence();

    void save_data_state(StateOut&);

    /// Set the current gradient and displacement.
    virtual void get_grad(const Ref<Function> &);
    virtual void get_x(const Ref<Function> &);
    virtual void set_nextx(const RefSCVector &);

    /// Set the current gradient and displacement to null.
    virtual void reset();

    /// Return nonzero if the optimization has converged.
    virtual int converged();

    void print(std::ostream& = ExEnv::out0()) const;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
