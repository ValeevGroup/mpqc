//
// point.cc
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
#pragma implementation
#endif

#include <stdlib.h>
#include <string.h>

#include <util/state/stateio.h>
#include <math/topology/point.h>
#include <math/scmat/matrix.h>
#include <util/keyval/keyval.h>
#include <util/misc/formio.h>

#define CLASSNAME Point
#define PARENTS public SavableState
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
Point::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

// Constructors
Point::Point(int in_dim) {dim=in_dim; x=new double[dim];}

Point::Point(const double *y, int in_dim)
{
    dim = in_dim;
    if (dim > 0) {
        x=new double[dim];
        if (!memcpy (x,y,dim*sizeof(double)))
            cerr << node0 << indent << "Bad copy in Point cnstrctr\n";
      }
}
    
Point::Point(const Point &in_point)
{
    dim=in_point.dim;
    if (dim > 0) {
        x=new double[dim];
        if (!memcpy (x,in_point.x,dim*sizeof(double)))
            cerr << node0 << indent << "Bad copy in Point copy ctor\n";
      }
}
    
Point::Point(RefSCVector &in)
{
    dim=in.dim().n();
    if (dim > 0) {
        x=new double[dim];
        in.convert(x);
      }
}

Point::Point(double a)
{
  dim = 1;
  x = new double[1];
  x[0] = a;
}

Point::Point(double a,double b)
{
  dim = 2;
  x = new double[2];
  x[0] = a;
  x[1] = b;
}

Point::Point(double a,double b,double c)
{
  dim = 3;
  x = new double[3];
  x[0] = a;
  x[1] = b;
  x[2] = c;
}

Point::Point(const RefKeyVal&k)
{
  dim = k->count();
  if (dim) x = new double[dim];
  for (int i=0; i<dim; i++) {
      x[i] = k->doublevalue(i);
    }
}

Point::Point(StateIn&s):
  SavableState(s)
{
  s.get(dim);
  if (dim > 0) s.get(x);
}

void
Point::save_data_state(StateOut& so)
{
  so.put(dim);
  if (dim > 0) so.put(x,dim);
}

// Destructor
Point::~Point(void) { if (dim > 0) delete[] x;};

// Copy to double *
double *Point::copy(void) const
{
  if (dim > 0) {
      double *x_out = new double[dim];
      if (!memcpy (x_out,x,dim*sizeof(double)))
          cerr << node0 << indent << "Bad copy in point copy function\n";
      return x_out;
    }
  return 0;
}

Point& Point::operator=(const Point&p)
{
  if (this != &p) {
      if (dim > 0) delete[] x;
      dim = p.dim;
      if (dim > 0) {
          x = new double[dim];
          for (int i=0; i<dim; i++) x[i] = p.x[i];
        }
    }
  return *this;
}

    // Set an element to a particular value
double &Point::operator[](int i)
{
    if (i < 0 || i >= dim)
    {
        cerr << node0 << indent << "Dimension out of bounds in Point[]\n";
	abort();
    }
    return x[i];
}

    // Set an element to a particular value
const double &Point::operator[](int i) const
{
    if (i < 0 || i >= dim)
    {
        cerr << node0 << indent << "Dimension out of bounds in Point[]\n";
	abort();
    }
    return x[i];
}

// Print out a Point
void Point::print(ostream& os)
{
    int i;
    for (i=0;i<dim;i++)	os << node0 << scprintf(" %12.8f",x[i]);
    os << node0 << endl;
}

void
Point::resize(int d)
{
  if (dim > 0) delete[] x;
  dim = d;
  if (dim > 0) x = new double[dim];
}

int Point::dimension() const
{
  return dim;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
