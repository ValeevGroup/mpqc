//
// pointset.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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

#ifndef _libQC_pointset_h
#define _libQC_pointset_h

#include <math/topology/point.h>
#include <util/container/ptrset.h>

template <class Type>
class PointSetElem<Type>
{
 private:
  Type obj;
  Point pnt;
 public:
  PointSetElem(Point p, Type t):pnt(p),obj(t) {};
  Point& point() { return pnt; };
  Type& object() { return obj; };
};

template <class Type>
class PointSet:
{
 private:
  PtrSet<PointSetElem<Type>> impl;
 public:
  PointSet() {};
  ~PointSet() {};
  length() { return impl.length(); }
  empty() { return 
};

#endif
