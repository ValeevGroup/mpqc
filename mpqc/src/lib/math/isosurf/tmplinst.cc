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
#include <scconfig.h>
#endif
#include <math/isosurf/shape.h>
#include <math/isosurf/triangle.h>

#ifdef HAVE_STL
#include <vector>
#endif

#ifdef EXPLICIT_TEMPLATE_INSTANTIATION

#ifdef HAVE_STL
template class vector<Ref<Vertex> >;
template class vector<Ref<Edge> >;
template class vector<Ref<Triangle> >;
#else
template class Array<Ref<Vertex> >;
template class Array<Ref<Edge> >;
template class Array<Ref<Triangle> >;
#endif

// Vertex
template class EAVLMMapNode<Ref<Vertex>, AVLMapNode<Ref<Vertex>, int> >;
template class EAVLMMap<Ref<Vertex>, AVLMapNode<Ref<Vertex>, int> >;
template class AVLMapNode<Ref<Vertex>, int>;
template class AVLMap<Ref<Vertex>, int>;
template class AVLSet<Ref<Vertex> >;

// Edge
template class EAVLMMapNode<Ref<Edge>, AVLMapNode<Ref<Edge>, int> >;
template class EAVLMMap<Ref<Edge>, AVLMapNode<Ref<Edge>, int> >;
template class AVLMapNode<Ref<Edge>, int>;
template class AVLMap<Ref<Edge>, int>;
template class AVLSet<Ref<Edge> >;

// Triangle
template class EAVLMMapNode<Ref<Triangle>, AVLMapNode<Ref<Triangle>, int> >;
template class EAVLMMap<Ref<Triangle>, AVLMapNode<Ref<Triangle>, int> >;
template class AVLMapNode<Ref<Triangle>, int>;
template class AVLMap<Ref<Triangle>, int>;
template class AVLSet<Ref<Triangle> >;

// Shape
template class EAVLMMapNode<Ref<Shape>, AVLMapNode<Ref<Shape>, int> >;
template class EAVLMMap<Ref<Shape>, AVLMapNode<Ref<Shape>, int> >;
template class AVLMapNode<Ref<Shape>, int>;
template class AVLMap<Ref<Shape>, int>;
template class AVLSet<Ref<Shape> >;

// (mixed)
template class EAVLMMapNode<Ref<Vertex>, AVLMapNode<Ref<Vertex>, AVLSet<Ref<Edge> > > >;
template class EAVLMMap<Ref<Vertex>, AVLMapNode<Ref<Vertex>, AVLSet<Ref<Edge> > > >;
template class AVLMapNode<Ref<Vertex>, AVLSet<Ref<Edge> > >;
template class AVLMap<Ref<Vertex>, AVLSet<Ref<Edge> > >;

template class EAVLMMapNode<Ref<Vertex>, AVLMapNode<Ref<Vertex>, Ref<Edge> > >;
template class EAVLMMap<Ref<Vertex>, AVLMapNode<Ref<Vertex>, Ref<Edge> > >;
template class AVLMapNode<Ref<Vertex>, Ref<Edge> >;
template class AVLMap<Ref<Vertex>, Ref<Edge> >;

template class EAVLMMapNode<Ref<Vertex>, AVLMapNode<Ref<Vertex>, AVLSet<Ref<Triangle> > > >;
template class EAVLMMap<Ref<Vertex>, AVLMapNode<Ref<Vertex>, AVLSet<Ref<Triangle> > > >;
template class AVLMapNode<Ref<Vertex>, AVLSet<Ref<Triangle> > >;
template class AVLMap<Ref<Vertex>, AVLSet<Ref<Triangle> > >;
#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
