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

#include <math/isosurf/shape.h>
#include <math/isosurf/triangle.h>

#ifdef __GNUG__
// Vertex
template class EAVLMMapNode<RefVertex, AVLMapNode<RefVertex, int> >;
template class EAVLMMap<RefVertex, AVLMapNode<RefVertex, int> >;
template class AVLMapNode<RefVertex, int>;
template class AVLMap<RefVertex, int>;
template class AVLSet<RefVertex>;

template class EAVLMMapNode<int, AVLMapNode<int, RefVertex> >;
template class EAVLMMap<int, AVLMapNode<int, RefVertex> >;
template class AVLMapNode<int, RefVertex>;
template class AVLMap<int, RefVertex>;

// Edge
template class EAVLMMapNode<RefEdge, AVLMapNode<RefEdge, int> >;
template class EAVLMMap<RefEdge, AVLMapNode<RefEdge, int> >;
template class AVLMapNode<RefEdge, int>;
template class AVLMap<RefEdge, int>;
template class AVLSet<RefEdge>;

template class EAVLMMapNode<int, AVLMapNode<int, RefEdge> >;
template class EAVLMMap<int, AVLMapNode<int, RefEdge> >;
template class AVLMapNode<int, RefEdge>;
template class AVLMap<int, RefEdge>;

// Triangle
template class EAVLMMapNode<RefTriangle, AVLMapNode<RefTriangle, int> >;
template class EAVLMMap<RefTriangle, AVLMapNode<RefTriangle, int> >;
template class AVLMapNode<RefTriangle, int>;
template class AVLMap<RefTriangle, int>;
template class AVLSet<RefTriangle>;

template class EAVLMMapNode<int, AVLMapNode<int, RefTriangle> >;
template class EAVLMMap<int, AVLMapNode<int, RefTriangle> >;
template class AVLMapNode<int, RefTriangle>;
template class AVLMap<int, RefTriangle>;

// Shape
template class EAVLMMapNode<RefShape, AVLMapNode<RefShape, int> >;
template class EAVLMMap<RefShape, AVLMapNode<RefShape, int> >;
template class AVLMapNode<RefShape, int>;
template class AVLMap<RefShape, int>;
template class AVLSet<RefShape>;

// (mixed)
template class EAVLMMapNode<RefVertex, AVLMapNode<RefVertex, AVLSet<RefEdge> > >;
template class EAVLMMap<RefVertex, AVLMapNode<RefVertex, AVLSet<RefEdge> > >;
template class AVLMapNode<RefVertex, AVLSet<RefEdge> >;
template class AVLMap<RefVertex, AVLSet<RefEdge> >;

template class EAVLMMapNode<RefVertex, AVLMapNode<RefVertex, RefEdge> >;
template class EAVLMMap<RefVertex, AVLMapNode<RefVertex, RefEdge> >;
template class AVLMapNode<RefVertex, RefEdge>;
template class AVLMap<RefVertex, RefEdge>;

template class EAVLMMapNode<RefVertex, AVLMapNode<RefVertex, AVLSet<RefTriangle> > >;
template class EAVLMMap<RefVertex, AVLMapNode<RefVertex, AVLSet<RefTriangle> > >;
template class AVLMapNode<RefVertex, AVLSet<RefTriangle> >;
template class AVLMap<RefVertex, AVLSet<RefTriangle> >;
#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
