//
// compute.h
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

#ifndef _util_misc_compute_h
#define _util_misc_compute_h

#include <stdio.h>
#include <util/container/avlset.h>
#include <util/state/state.h>

class ResultInfo;
class StateIn;
class StateOut;

typedef ResultInfo* ResultInfoP;

class Compute
{
   friend class ResultInfo;
   friend class AccResultInfo;
  private:
    AVLSet<ResultInfoP> _results;
    void add(ResultInfo*);

    // Prohibit copy
    Compute(const Compute&) {};

  protected:
    // Recompute at least the results that have compute true
    // and are not already computed.  This should only be called
    // by Result's members.
    virtual void compute() = 0;
  public:
    Compute();
    virtual ~Compute();

    // Marks all results as being out of date.  Any subsequent access
    // to results will cause Compute::compute() to be called.
    virtual void obsolete();
};

// Usually Result_dec(Type) will be used to create a result that has
// a particular datum associated with it, however simple Result's can
// also be declared to keep track of datum's for which it is awkward
// to use Result_dec.
class ResultInfo
{
  protected:
    int _compute;
    int _computed;
    Compute* _c;
    // This make sure that the datum is up to date.  If it is not then
    // Compute::compute() will be called.
    virtual void update();
  protected:
    ResultInfo(StateIn&,Compute*);
    ResultInfo(const ResultInfo&,Compute*);
    virtual void save_data_state(StateOut&);
    virtual void restore_state(StateIn&);
    ResultInfo& operator=(const ResultInfo&);
  public:
    ResultInfo(Compute*c);
    virtual ~ResultInfo();
    int& compute() { return _compute; }
    int compute(int c) { int r = _compute; _compute = c; return r; }
    int& computed() { return _computed; }
    virtual int needed();
};

// This is like result but the accuracy with which a result was computed as
// well as the desired accuracy are stored.  A _computed datum may has an
// actual accuracy greater than, equal to, or less than the computed
// accuracy.
class AccResultInfo: public ResultInfo
{
  private:
    double _actual_accuracy;
    double _desired_accuracy;
  protected:
    AccResultInfo(StateIn&,Compute*);
    AccResultInfo(const AccResultInfo&,Compute*);
    virtual void save_data_state(StateOut&);
    virtual void restore_state(StateIn&);
    AccResultInfo& operator=(const AccResultInfo&);
    void update();
  public:
    AccResultInfo(Compute*c);
    ~AccResultInfo();
    double actual_accuracy() const;
    double desired_accuracy() const;
    void set_desired_accuracy(double);
    void set_actual_accuracy(double);
    int computed_to_desired_accuracy()
        { return computed() && _actual_accuracy <= _desired_accuracy; }
    int needed();
};

#include <util/misc/comptmpl.h>

typedef NCResult<int> Resultint;
typedef NCResult<double> Resultdouble;
typedef NCAccResult<double> AccResultdouble;

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
