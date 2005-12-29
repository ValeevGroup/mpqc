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

namespace sc {

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
    int n() { return ((l_>=0)?((((l_)+2)*((l_)+1))>>1):0); }
    /// Returns the exponent of x.
    int a() { return a_; }
    /// Returns the exponent of y.
    int b() { return b_; }
    /// Returns the exponent of z.
    int c() { return c_; }
    /// Returns the angular momentum.
    int l() { return l_; }
    /// Returns a() if i==0, b() if i==1, and c() if i==2.
    int l(int i) { return i ? (i==1 ? b_ : c_) : a_; }
    /** Returns the number of the current basis function within the shell.
        This starts at 0 and sequentially increases as next() is called. */
    int bfn() { return bfn_; }
};

/** RedundantCartesianIter objects loop through all possible combinations
    of a given number of axes.  This is used to compute the transformation
    matrices that maps a set of Cartesian functions to another set of
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
    operator int() { return !done_; }

    /// The current exponent of x.
    int a();
    /// The current exponent of y.
    int b();
    /// The current exponent of z.
    int c();
    /// The angular momentum.
    int l() { return l_; }
    /// Returns a() if i==0, b() if i==1, and c() if i==2.
    int l(int i);
    /// Return the i'th axis.
    int axis(int i) { return axis_[i]; }
};

inline void
RedundantCartesianIter::start()
{
  if (l_==0)
    done_ = 1;
  else
    done_ = 0;

  for (int i=0; i<l_; i++)
    axis_[i] = 0;
}

inline void
RedundantCartesianIter::next()
{
  for (int i=0; i<l_; i++) {
    if (axis_[i] == 2)
      axis_[i] = 0;
    else {
      axis_[i]++;
      return;
    }
  }
  done_ = 1;
}

inline int
RedundantCartesianIter::l(int axis)
{
  int i;
  int r = 0;
  for (i=0; i<l_; i++) if (axis_[i]==axis) r++;
  return r;
}

inline int
RedundantCartesianIter::a()
{
  return l(0);
}

inline int
RedundantCartesianIter::b()
{
  return l(1);
}

inline int
RedundantCartesianIter::c()
{
  return l(2);
}

/** Like RedundantCartesianIter, except a, b, and c are fixed to a given
    value. */
class RedundantCartesianSubIter {
  private:
    int done_;
    int l_;
    int e_[3];
    int *axis_;

    // the locations of the z's in the axis array
    int *zloc_;
    // the locations of the y's in the subarray after the z's are removed
    int *yloc_;

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
    operator int() const { return !done_; }

    /// The current exponent of x.
    int a() const { return e_[0]; }
    /// The current exponent of y.
    int b() const { return e_[1]; }
    /// The current exponent of z.
    int c() const { return e_[2]; }
    /// The angular momentum.
    int l() const { return l_; }
    /// Returns a() if i==0, b() if i==1, and c() if i==2.
    int l(int i) { return e_[i]; }
    /// Return the i'th axis.
    int axis(int i) { return axis_[i]; }
};

}
  
#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
