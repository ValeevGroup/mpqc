//
// vertex.cc
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

#include <util/misc/formio.h>
#include <math/isosurf/vertex.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////////
// Vertex (a point and a gradient)

Vertex::Vertex(const SCVector3&point):
  _point(point),
  _normal(0)
{
}

Vertex::Vertex(const SCVector3&point,const SCVector3&normal):
  _point(point),
  _normal(new SCVector3(normal))
{
  double dot = _normal->dot(*_normal);
  if (dot < 0.999999 || dot > 1.000001) {
      ExEnv::outn() << "Vertex: ctor: bad normal\n" << endl;
      abort();
    }
}

Vertex::Vertex():
  _normal(0)
{
}

Vertex::~Vertex()
{
  if (_normal) delete _normal;
}

void
Vertex::set_point(const SCVector3&p)
{
  _point = p;
}

void
Vertex::set_normal(const SCVector3&p)
{
  if (_normal) {
      *_normal = p;
    }
  else {
      _normal = new SCVector3(p);
    }
  double dot = _normal->dot(*_normal);
  if (dot < 0.999999 || dot > 1.000001) {
      ExEnv::outn() << "Vertex::set_normal: bad normal\n" << endl;
      abort();
    }
}

Vertex::operator SCVector3&()
{
  return _point;
}

void
Vertex::print(ostream&o)
{
  int i;
  o << indent << "Vertex:";
  for (i=0; i<3; i++)  {
      o << scprintf(" %8.5f", _point[i]);
    }
  if (_normal) {
      for (i=0; i<3; i++)  {
          o << scprintf(" %8.5f", normal()[i]);
        }
    }
  o << endl;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
