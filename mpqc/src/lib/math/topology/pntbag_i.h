//
// pntbag_i.h
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

#ifdef __GNUC__
#pragma interface
#endif

#ifdef INLINE_FUNCTIONS
#define INLINE inline
#else
#define INLINE
#endif

INLINE int PointBag_double::length()
{
  return impl.length();
}

INLINE void PointBag_double::add(const Point& p, double x)
{
  impl.add((void*)new PointBagElem_double(p,x));
}

INLINE void PointBag_double::add(const PointBagElem_double& e)
{
  impl.add((void*)new PointBagElem_double(e));
}

INLINE int PointBag_double::owns(Pix i)
{
  return impl.owns(i);
}

INLINE double& PointBag_double::operator()(Pix i)
{
  return ((PointBagElem_double*)impl(i).getptr())->object();
}

INLINE double& PointBag_double::get(Pix i)
{
  return operator()(i);
}

INLINE Point& PointBag_double::point(Pix i)
{
  return ((PointBagElem_double*)impl(i).getptr())->point();
}

INLINE Pix PointBag_double::first()
{
  return impl.first();
}

INLINE void PointBag_double::next(Pix&i)
{
  impl.next(i);
}

#undef INLINE
