//
// function.h
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

#ifdef __GNUC__
#pragma interface
#endif

#ifndef _math_optimize_function_h
#define _math_optimize_function_h

#include <math.h>
#include <float.h>

#include <util/state/state.h>
#include <math/optimize/transform.h>
#include <math/scmat/matrix.h>
#include <math/scmat/result.h>

//. The \clsnm{Function} class is an abstract base class that,
//given a set of coordinates, will compute a value and possibly
//a gradient and hessian at that point.
class Function: virtual_base public SavableState, public Compute {
#   define CLASSNAME Function
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    RefSCMatrixKit matrixkit_;

    RefSCVector x_;                     // variables
    RefSCDimension dim_;                // dimension of x_
    AccResultdouble value_;             // value of function at x_
    AccResultRefSCVector gradient_;     // gradient at x_
    AccResultRefSymmSCMatrix hessian_;  // hessian at x_

    //. Update the various computable results.
    virtual void set_value(double);
    virtual void set_gradient(RefSCVector&);
    virtual void set_hessian(RefSymmSCMatrix&);

    //. Set the \clsnmref{SCMatrixKit} that should be used to
    //construct the requisite vectors and matrices.
    virtual void set_matrixkit(const RefSCMatrixKit&);
    virtual void set_dimension(const RefSCDimension&);

    //. Set the accuracies with which the various computables
    //must be computed.
    virtual void set_actual_value_accuracy(double);
    virtual void set_actual_gradient_accuracy(double);
    virtual void set_actual_hessian_accuracy(double);

    //. Get read/write access to the coordinates for modification.
    RefSCVector& get_x_reference() { obsolete(); return x_; }

    void do_change_coordinates(const RefNonlinearTransform&);
  public:
    //. The standard constructors and destructor.
    Function();
    Function(StateIn&);
    Function(const Function&);
    Function(const RefKeyVal&, double funcacc = DBL_EPSILON,
             double gradacc = DBL_EPSILON, double hessacc = DBL_EPSILON);
    virtual ~Function();

    Function & operator=(const Function&);

    //. Return the \clsnmref{SCMatrixKit} used to construct
    //vectors and matrices.
    RefSCMatrixKit matrixkit() const;
    //. Return the \clsnmref{SCDimension} of the problem.
    RefSCDimension dimension() const;

    virtual void save_data_state(StateOut&);

    //. Return the value of the function.
    virtual double value();
    //. Returns nonzero if the current value is not up-to-date.
    int value_needed() const;
    //. If passed a nonzero number, compute the value the next
    //time \srccd{compute()} is called.  Return a nonzero number
    //if the value was previously to be computed.
    int do_value(int);
    AccResultdouble& value_result() { return value_; }

    //. Set the accuracy to which the value is to be computed.
    virtual void set_desired_value_accuracy(double);
    //. Return the accuracy with which the value has been computed.
    virtual double actual_value_accuracy() const;
    //. Return the accuracy with which the value is to be computed.
    virtual double desired_value_accuracy() const;

    //. These are analogous to the routines that deal with values,
    //but work with gradients instead.
    virtual RefSCVector gradient();
    int gradient_needed() const;
    int do_gradient(int);
    virtual void set_desired_gradient_accuracy(double);
    virtual double actual_gradient_accuracy() const;
    virtual double desired_gradient_accuracy() const;
    AccResultRefSCVector& gradient_result() { return gradient_; }

    //. These are analogous to the routines that deal with values,
    //but work with the hessian instead.
    virtual RefSymmSCMatrix hessian();
    int hessian_needed() const;
    int do_hessian(int);
    virtual void set_desired_hessian_accuracy(double);
    virtual double actual_hessian_accuracy() const;
    virtual double desired_hessian_accuracy() const;
    AccResultRefSymmSCMatrix& hessian_result() { return hessian_; }

    // hessian by gradients at finite displacements
    // virtual RefSCMatrix fd1_hessian();

    //. Compute a quick, approximate hessian.
    virtual void guess_hessian(RefSymmSCMatrix&);
    virtual RefSymmSCMatrix inverse_hessian(RefSymmSCMatrix&);

    //. Information about the availability of values, gradients,
    //and hessians.
    virtual int value_implemented() const;
    virtual int gradient_implemented() const;
    virtual int hessian_implemented() const;

    //. Set and retrieve the coordinate values.
    virtual void set_x(const RefSCVector&);
    RefSCVector get_x() const { return x_.copy(); }
    const RefSCVector& get_x_no_copy() const { return x_; }

    //. An optimizer can call change coordinates periodically
    //to give the function an opportunity to change its coordinate
    //system.  A return value of 0 means the coordinates were not
    //changed.  Otherwise, a transform object to the new coordinate
    //system is return.  The function object applies the transform
    //to any objects it contains.  This will obsolete the function
    //data.
    virtual RefNonlinearTransform change_coordinates();

    //. Print information about the object.
    virtual void print(ostream& = cout) const;
};
SavableState_REF_dec(Function);

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
