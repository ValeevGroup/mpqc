//
// pointbag.cc
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

#include <math/topology/pointbag.h>

PointBagElem_double::PointBagElem_double(const Point&p, const double&t):
  obj(t),
  pnt(p)
{
}
PointBagElem_double::PointBagElem_double(const PointBagElem_double&b):
  obj(b.obj),
  pnt(b.pnt)
{
}
PointBagElem_double::~PointBagElem_double()
{
}

Point& PointBagElem_double::point()
{
  return pnt;
}

double& PointBagElem_double::object()
{
  return obj;
}

///////////////////////////////////////////////////////////////////

PointBag_double::PointBag_double()
{
}
PointBag_double::PointBag_double(PointBag_double&b)
{
  for (Pix i=b.first(); i!=0; b.next(i)) {
      add(*((PointBagElem_double*) &b.impl(i)));
    }
}
PointBag_double::~PointBag_double()
{
  for (Pix i=first(); i!=0; next(i)) {
      delete (PointBagElem_double*) impl(i).getptr();
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
