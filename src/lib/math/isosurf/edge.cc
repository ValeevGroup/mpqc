//
// edge.cc
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

#include <math/isosurf/triangle.h>
#include <math/isosurf/edge.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////////
// Edge

Edge::Edge(const Ref<Vertex> &p1,
     const Ref<Vertex> &p2,
     const Ref<Vertex> &p3)
{
  _order = 2;
  _vertices = new Ref<Vertex>[3];
  _vertices[0]=p1; _vertices[1]=p2; _vertices[2]=p3;
}

Edge::Edge(const Ref<Vertex> &p1,
     const Ref<Vertex> &p2,
     const Ref<Vertex> &p3,
     const Ref<Vertex> &p4)
{
  _order = 3;
  _vertices = new Ref<Vertex>[4];
  _vertices[0]=p1; _vertices[1]=p2; _vertices[2]=p3; _vertices[3]=p4;
}

Edge::~Edge()
{
  delete[] _vertices;
}

double Edge::straight_length()
{
  SCVector3 BA = vertex(1)->point() - vertex(0)->point();
  return BA.norm();
}

void Edge::add_vertices(std::set<Ref<Vertex> >&set)
{
  set.insert(_vertices[0]);
  set.insert(_vertices[_order]);
}

void
Edge::set_order(int order, const Ref<Volume>&vol,double isovalue)
{
  Ref<Vertex> *newvertices = new Ref<Vertex>[order+1];
  newvertices[0] = vertex(0);
  newvertices[order] = vertex(1);
  delete[] _vertices;
  _vertices = newvertices;
  _order = order;

  SCVector3 pv[2];
  SCVector3 norm[2];

  int i;
  for (i=0; i<2; i++) {
      pv[i] = vertex(i)->point();
      norm[i] = vertex(i)->normal();
    }

  for (i=1; i<_order; i++) {
      double I = (1.0*i)/order;
      double J = (1.0*(_order - i))/order;
      SCVector3 interpv = I * pv[0] + J * pv[1];
      SCVector3 start(interpv);
      SCVector3 interpnorm = I * norm[0] + J * norm[1];
      SCVector3 newpoint;
      vol->solve(start,interpnorm,isovalue,newpoint);
      vol->set_x(newpoint);
      if (vol->gradient_implemented()) {
          vol->get_gradient(interpnorm);
        }
      interpnorm.normalize();
      _vertices[i] = new Vertex(newpoint,interpnorm);
    }

}

int
Edge::interpolate(double r, SCVector3&point, SCVector3&norm)
{
  int i;
  double s = 1.0 - r;
  double rcoef[Triangle::max_order+1];
  double scoef[Triangle::max_order+1];
  double spacing = 1.0/_order;
  rcoef[0] = 1.0;
  scoef[0] = 1.0;
  for (i=1; i<=_order; i++) {
      rcoef[i] = rcoef[i-1]*(r-(i-1)*spacing)/(i*spacing);
      scoef[i] = scoef[i-1]*(s-(i-1)*spacing)/(i*spacing);
    }
  int has_norm = 1;
  for (i=0; i<=_order; i++) {
      if ((has_norm = (has_norm && _vertices[i]->has_normal()))) break;
    }
  point = 0.0;
  norm = 0.0;
  for (i=0; i<=_order; i++) {
      int j=_order-i;
      double coef = rcoef[i]*scoef[j];
      point += coef * _vertices[i]->point();
      if (has_norm) norm += coef * _vertices[i]->normal();
    }
  if (has_norm) norm.normalize();
  return has_norm;
}

int
Edge::interpolate(double r, SCVector3&point, SCVector3&norm,
                  const Ref<Volume> &vol, double isovalue)
{
  // first guess
  int has_norm = interpolate(r,point,norm);

  if (!has_norm) {
      ExEnv::errn() << "Edge::interpolate with volume requires norm"
           << endl;
      abort();
    }

  // refine guess
  SCVector3 newpoint;
  vol->solve(point,norm,isovalue,newpoint);
  // compute the true normal
  vol->set_x(newpoint);
  if (vol->gradient_implemented()) {
      vol->get_gradient(norm);
      norm.normalize();
    }

  return has_norm || vol->gradient_implemented();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
