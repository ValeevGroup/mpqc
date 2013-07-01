//
// tmplinst.cc - template instantations for isosurf
//
// Copyright (C) 1998 Limit Point Systems, Inc.
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

#ifdef HAVE_CONFIG_H
#include <mpqc_config.h>
#endif
#include <math/isosurf/shape.h>
#include <math/isosurf/triangle.h>

#include <vector>

using namespace sc;

#ifdef EXPLICIT_TEMPLATE_INSTANTIATION

template class vector<Ref<Vertex> >;
template class vector<Ref<Edge> >;
template class vector<Ref<Triangle> >;

// Vertex
template class map<Ref<Vertex>, int>;
template class set<Ref<Vertex> >;

// Edge
template class map<Ref<Edge>, int>;
template class set<Ref<Edge> >;

// Triangle
template class map<Ref<Triangle>, int>;
template class set<Ref<Triangle> >;

// Shape
template class map<Ref<Shape>, int>;
template class set<Ref<Shape> >;

// (mixed)
template class map<Ref<Vertex>, set<Ref<Edge> > >;
template class map<Ref<Vertex>, Ref<Edge> >;
template class map<Ref<Vertex>, set<Ref<Triangle> > >;
#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
