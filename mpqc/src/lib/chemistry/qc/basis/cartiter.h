//
// cartiter.h
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

#ifndef _chemistry_qc_basis_cartiter_h
#define _chemistry_qc_basis_cartiter_h

#ifdef __GNUC__
#pragma interface
#endif

/** CartesianIter gives the ordering of the Cartesian functions
    within a shell for the particular integrals specialization. */
class CartesianIter {
  protected:
    int a_;
    int b_;
    int c_;
    int l_;
    int bfn_;

  public:
    /// Initialize an iterator for the given angular momentum.
    CartesianIter(int l);
    virtual ~CartesianIter();

    /// Start the iteration.
    virtual void start() =0;
    /// Move to the next Cartesian function.
    virtual void next() =0;
    /// Returns nonzero if the iterator currently hold valid data.
    virtual operator int() =0;

    /// Returns the number of Cartesian functions.
    int n();
    /// Returns the exponent of x.
    int a();
    /// Returns the exponent of y.
    int b();
    /// Returns the exponent of z.
    int c();
    /// Returns the angular momentum.
    int l();
    /// Returns a() if i==0, b() if i==1, and c() if i==2.
    int l(int i);
    /** Returns the number of the current basis function within the shell.
        This starts at 0 and sequentially increases as next() is called. */
    int bfn();
};

/** RedundantCartesianIter objects loop through all possible combinations
    of a given number of axes.  This is used to compute the transformation
    matrices that maps a set of Cartesian functions into to another set of
    Cartesian functions in a rotated coordinate system. */
class RedundantCartesianIter {
  private:
    int done_;
    int l_;
    int *axis_;

  public:
    /// Create a object for the given angular momentum.
    RedundantCartesianIter(int l);
    virtual ~RedundantCartesianIter();

    /// Return the current Cartesian basis function number.
    virtual int bfn() =0;

    /// Initialize the iterator.
    void start();
    /// Move to the next combination of axes.
    void next();
    /// Returns nonzero if the iterator currently hold valid data.
    operator int();

    /// The current exponent of x.
    int a();
    /// The current exponent of y.
    int b();
    /// The current exponent of z.
    int c();
    /// The angular momentum.
    int l();
    /// Returns a() if i==0, b() if i==1, and c() if i==2.
    int l(int i);
    /// Return the i'th axis.
    int axis(int i);
};

/** Like RedundantCartesianIter, except a, b, and c are fixed to a given
    value. */
class RedundantCartesianSubIter {
  private:
    int done_;
    int l_;
    int e_[3];
    int *axis_;

    void advance();
    int valid();

  public:
    /// Create a object for the given angular momentum.
    RedundantCartesianSubIter(int l);
    virtual ~RedundantCartesianSubIter();

    /// Return the current Cartesian basis function number.
    virtual int bfn() =0;

    /** Initialize the iterator.  The constraints on a, b, and c are
        given as arguments. */
    void start(int a, int b, int c);
    /// Move to the next combination of axes.
    void next();
    /// Returns nonzero if the iterator currently hold valid data.
    operator int();

    /// The current exponent of x.
    int l();
    /// The current exponent of y.
    int a();
    /// The current exponent of z.
    int b();
    /// The angular momentum.
    int c();
    /// Returns a() if i==0, b() if i==1, and c() if i==2.
    int l(int i);
    /// Return the i'th axis.
    int axis(int i);
};

#ifdef INLINE_FUNCTIONS
#include <chemistry/qc/basis/cartiter_i.h>
#endif
  
#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
