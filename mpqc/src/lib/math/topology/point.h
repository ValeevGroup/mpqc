//
// point.h
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

#ifndef _libQC_point_h
#define _libQC_point_h

#ifdef __GNUC__
#pragma interface
#endif

// This class implements arbitrary dimensional point using double precision #'s
#include <iostream.h>
#include <string.h>
#include <ctype.h>

#include <util/class/class.h>
#include <util/state/state.h>
#include <util/container/set.h>
#include <util/container/array.h>

class RefSCVector;
class RefKeyVal;

class Point : public SavableState
{
#   define CLASSNAME Point
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    int dim;
    double *x;
  
  public:
  // Constructors
    Point(int in_dim=3);
    Point(const double *y, int in_dim=3);
    Point(const Point &in_point);
    Point(RefSCVector&);
    Point(double);
    Point(double,double);
    Point(double,double,double);
    Point(StateIn&);
    Point(const RefKeyVal&);
  
  // Destructor
    ~Point(void);
  
    int dimension() const;
    void resize(int dim);

  // Copy to double *
    double *copy(void) const;

    Point& operator=(const Point&);
  
  // Set an element to a particular value
    double &operator[](int i);
    const double &operator[](int i) const;
  
  // Print out a Point
    void print(ostream& = cout);

  // save and restore state
    void save_data_state(StateOut&);
    void restore_data_state(int,StateIn&);

};

DescribedClass_REF_dec(Point);
ARRAY_dec(RefPoint);
SET_dec(RefPoint);
ARRAYSET_dec(RefPoint);

class cart_point {
 private:
  double r[3];
 public:
  cart_point();
  cart_point(const cart_point&p);
  cart_point(const double*);
  ~cart_point();
  double& operator[](int i);
  const double& operator[](int i) const;
  double& x();
  double& y();
  double& z();
  const double& x() const;
  const double& y() const;
  const double& z() const;
};

typedef cart_point Point3;

typedef struct {
  double r;
  double theta;
  double phi;
} sph_point;

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
