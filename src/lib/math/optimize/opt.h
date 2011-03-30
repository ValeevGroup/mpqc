//
// opt.h
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

#ifndef _math_optimize_opt_h
#define _math_optimize_opt_h

#include <util/group/message.h>
#include <util/state/state.h>
#include <util/class/class.h>
#include <math/scmat/matrix.h>
#include <math/optimize/function.h>
#include <math/optimize/conv.h>

namespace sc {

// //////////////////////////////////////////////////////////////////////

/** The Optimize class is an abstract base class for classes
    that find the extreme points of Function's. */
class Optimize: virtual public SavableState {
  protected:
    int max_iterations_;
    int n_iterations_;
    int ckpt_;
    int print_timings_;
    double max_stepsize_;
    std::string ckpt_file_;
    Ref<Function> function_;
    Ref<Convergence> conv_;
    Ref<MessageGrp> msg_;
  public:
    Optimize();
    /// Restore the state of a Function object.
    Optimize(StateIn&);

    /** The KeyVal constructor reads the following information:

        <dl>

        <dt><tt>checkpoint</tt><dd> If true, the optimization will be
        checkpointed.  The default is false.

        <dt><tt>checkpoint_file</tt><dd> The name of the checkpoint file.
        The name defaults to opt_ckpt.dat.

        <dt><tt>max_iterations</tt><dd> The maximum number of interations.
        The default is 10.

        <dt><tt>max_stepsize</tt><dd> The maximum stepsize.  The default is
        0.6.

        <dt><tt>function</tt><dd> A Function object.  There is no default.

        <dt><tt>convergence</tt><dd> This can be either a floating point
        number or a Convergence object.  If it is a floating point number
        then it is the convergence criterion.  See the description
        Convergence class for the default.

        </dl> */
    Optimize(const Ref<KeyVal>&);
    virtual ~Optimize();

    void save_data_state(StateOut&);

    /** Do the optimization.  Returns nonzero if the optimization
        is complete. */
    virtual int optimize();

    /// Set up for checkpointing.
    void set_checkpoint();
    void set_checkpoint_file(const char*);

    /// Set the function to be optimized
    void set_function(const Ref<Function>&);
    
    /// Set the iteration limit.
    void set_max_iterations(int);
  
    /// Initialize the optimizer.
    virtual void init();
    /** Take a step.  Returns 1 if the optimization has converged,
        otherwise 0. */
    virtual int update() = 0;

    virtual void apply_transform(const Ref<NonlinearTransform>&);

    /// Returns information about the Function being optimized.
    Ref<Function> function() const { return function_; }
    Ref<SCMatrixKit> matrixkit() const { return function_->matrixkit(); }
    RefSCDimension dimension() const { return function_->dimension(); }

    void print(std::ostream& = ExEnv::out0()) const;
};


/** The LineOpt abstract class is used to perform one dimensional
optimizations.*/
class LineOpt: public Optimize {

  protected:

    RefSCVector initial_x_;
    double initial_value_;
    RefSCVector initial_grad_;
    RefSCVector search_direction_;
    Ref<Function> function_;

  public:

    LineOpt();
    LineOpt(StateIn&);
    LineOpt(const Ref<KeyVal>&);
    ~LineOpt();
    void save_data_state(StateOut&);

    /** Initializes the line search object. Argument is a search direction.
      * Use of this method assumes the Optimize base class already has a 
      * function object (got it from a keyval or elsewhere). */
    virtual void init(RefSCVector& direction);
    /** Initializes the line search object. First argument is a search 
      * direction, second argument is a function object to optimize.
      * Use this method when a function must be passed to the Optimize 
      * base class. */
    virtual void init(RefSCVector& direction, Ref<Function> function);
    /// Applies a nonlinear transform.
    void apply_transform(const Ref<NonlinearTransform>&);

    void print(std::ostream& = ExEnv::out0()) const;
};

class Backtrack: public LineOpt {

 protected:
   double decrease_factor_;
   double backtrack_factor_;
   int force_search_;
    
   int sufficient_decrease(RefSCVector& step);

 public:
   Backtrack();
   Backtrack(const Ref<KeyVal>&);
   Backtrack(StateIn&s);
   ~Backtrack();
   int update();
   void save_data_state(StateOut&s);
  
   int force_search() const { return force_search_; }
   /// Returns factor for sufficient decrease test
   double decrease_factor() { return decrease_factor_; }
   /// Sets factor for sufficient decrease test
   double set_decrease_factor( double factor ) 
   { double temp = decrease_factor_; decrease_factor_ = factor; return temp; }

   void print(std::ostream& = ExEnv::out0()) const;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
