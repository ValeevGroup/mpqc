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

#ifdef __GNUC__
#pragma interface
#endif

#include <util/class/class.h>
#include <util/state/state.h>
#include <util/keyval/keyval.h>

SavableState_REF_fwddec(SCExtrapData);
class SCExtrapData: public SavableState {
#   define CLASSNAME SCExtrapData
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SCExtrapData();
    SCExtrapData(StateIn&);
    virtual ~SCExtrapData();

    void save_data_state(StateOut&);
    
    virtual SCExtrapData* copy() = 0;
    virtual void zero() = 0;
    virtual void accumulate_scaled(double scale, const RefSCExtrapData&) = 0;
};
SavableState_REF_dec(SCExtrapData);

SavableState_REF_fwddec(SCExtrapError);
class SCExtrapError: public SavableState {
#   define CLASSNAME SCExtrapError
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SCExtrapError();
    SCExtrapError(StateIn&);
    virtual ~SCExtrapError();

    void save_data_state(StateOut&);
    
    virtual double error() = 0;
    virtual double scalar_product(const RefSCExtrapError&) = 0;
};
SavableState_REF_dec(SCExtrapError);

class SelfConsistentExtrapolation: public SavableState {
#   define CLASSNAME SelfConsistentExtrapolation
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  private:
    double error_;
    int errorset_;
    double tolerance_;
  protected:
    void set_error(double e) { error_ = e; errorset_ = 1; }
  public:
    SelfConsistentExtrapolation();
    SelfConsistentExtrapolation(StateIn&);
    SelfConsistentExtrapolation(const RefKeyVal&);
    ~SelfConsistentExtrapolation();

    void save_data_state(StateOut&);
    
    void set_tolerance(double t) { tolerance_ = t; }
    double tolerance() { return tolerance_; }
    double error() { return error_; }

    int converged() { return errorset_? error_ <= tolerance_ : 0; }

    // Makes a copy of data and returns the extrapolation in
    // data.  A reference to error is saved so a copy must
    // be given to extrapolate if error could be changed.
    virtual int extrapolate(const RefSCExtrapData& data,
                            const RefSCExtrapError& error) = 0;

    // Extrapolation should be started when this is called,
    // if it hasn't already started.  The default starting
    // point is implemenation dependent.  This member might
    // do nothing in some implementations.
    virtual void start_extrapolation();

    virtual void reinitialize() =0;
};
SavableState_REF_dec(SelfConsistentExtrapolation);

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
