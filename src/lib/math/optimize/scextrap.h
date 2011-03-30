//
// scextrap.h
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

#ifndef _math_optimize_scextrap_h
#define _math_optimize_scextrap_h

#include <util/class/class.h>
#include <util/state/state.h>
#include <util/keyval/keyval.h>

namespace sc {

/** SCExtrapData hold the data to be extrapolated needed by
    SelfConsistentExtrapolation.  */
class SCExtrapData: public SavableState {
  public:
    /// Construct a new SCExtrapData.
    SCExtrapData();
    /// Constructor to restore SCExtrapData from a StateIn object.
    SCExtrapData(StateIn&);
    virtual ~SCExtrapData();

    void save_data_state(StateOut&);
    
    /** Return a copy of this. */
    virtual SCExtrapData* copy() = 0;
    /** Set this to zero. */
    virtual void zero() = 0;
    /** The passed SCExtrapData is accumulated into this after applying the
        scaling factor. */
    virtual void accumulate_scaled(double scale, const Ref<SCExtrapData>&) = 0;
};


/** SCExtrapError holds the error data needed by SelfConsistentExtrapolation.
 */
class SCExtrapError: public SavableState {
  public:
    /// Construct a new SCExtrapError.
    SCExtrapError();
    /// Constructor to restore SCExtrapError from a StateIn object.
    SCExtrapError(StateIn&);
    virtual ~SCExtrapError();

    void save_data_state(StateOut&);
    
    /// Returns some measure of the total error.
    virtual double error() = 0;
    /// Performs a scalar product between this and the given SCExtrapError.
    virtual double scalar_product(const Ref<SCExtrapError>&) = 0;
};


/** The SelfConsistentExtrapolation abstract class is used to iteratively
solve equations requiring a self consistent solution, such as,

\f[ \bar{x}' = f(\bar{x}) \f]
*/
class SelfConsistentExtrapolation: public SavableState {
  private:
    double error_;
    int errorset_;
    double tolerance_;
  protected:
    void set_error(double e) { error_ = e; errorset_ = 1; }
  public:
    SelfConsistentExtrapolation();
    SelfConsistentExtrapolation(StateIn&);
    /** The only keyword read is <tt>tolerance</tt>, which is usually not needed
        since the objects using SelfConsistentExtrapolation should set the
        tolerances as needed.  */
    SelfConsistentExtrapolation(const Ref<KeyVal>&);
    ~SelfConsistentExtrapolation();

    void save_data_state(StateOut&);
    
    void set_tolerance(double t) { tolerance_ = t; }
    double tolerance() { return tolerance_; }
    double error() { return error_; }

    int converged() { return errorset_? error_ <= tolerance_ : 0; }

    // Makes a copy of data and returns the extrapolation in
    // data.  A reference to error is saved so a copy must
    // be given to extrapolate if error could be changed.
    virtual int extrapolate(const Ref<SCExtrapData>& data,
                            const Ref<SCExtrapError>& error) = 0;

    // Extrapolation should be started when this is called,
    // if it hasn't already started.  The default starting
    // point is implemenation dependent.  This member might
    // do nothing in some implementations.
    virtual void start_extrapolation();

    /// This must be called if the extrapolation object is to be used
    /// again. It should also be called before the first use of the object,
    /// if initial data needs to be given to the algorithm.  The data
    /// object will be copied.
    virtual void reinitialize(Ref<SCExtrapData> data=0) =0;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
