
extern "C" {
#include <stdio.h>
#include <math.h>
#include "p3dgen.h"
}
#include "surf.h"
#include "isosurf.h"

/////////////////////////////////////////////////////////////////////////
// Vertex (a point and a gradient)

REF_def(Vertex);
ARRAY_def(RefVertex);
SET_def(RefVertex);
ARRAYSET_def(RefVertex);
ARRAY_def(ArraysetRefVertex);

Vertex::Vertex(RefPoint&point,DVector&gradient):
  _point(point),
  _gradient(gradient)
{
}

Vertex::Vertex()
{
}

Vertex::~Vertex()
{
}

void
Vertex::print(FILE*fp)
{
  int i;
  fprintf(fp, "Vertex:");
  for (i=0; i<_point->dimension(); i++)  {
      fprintf(fp," %8.5f",_point->operator[](i));
    }
  for (i=0; i<_gradient.dim(); i++)  {
      fprintf(fp," %8.5f",_gradient[i]);
    }
  fprintf(fp,"\n");
}

/////////////////////////////////////////////////////////////////////////
// Edge

REF_def(Edge);
ARRAY_def(RefEdge);
SET_def(RefEdge);
ARRAYSET_def(RefEdge);
ARRAY_def(ArraysetRefEdge);

Edge::~Edge()
{
}

double Edge::length()
{
  DVector A(_vertices[0]->point());
  DVector B(_vertices[1]->point());
  DVector BA = B - A;
  return sqrt(BA.dot(BA));
}

void Edge::add_vertices(SetRefVertex&set)
{
  set.add(_vertices[0]);
  set.add(_vertices[1]);
}

/////////////////////////////////////////////////////////////////////////
// Edge4

REF_def(Edge4);
ARRAY_def(RefEdge4);
SET_def(RefEdge4);
ARRAYSET_def(RefEdge4);
ARRAY_def(ArraysetRefEdge4);

Edge4::~Edge4()
{
}

Edge4::Edge4(RefVertex p1,RefVertex p2,RefVolume vol,double isovalue):
  Edge(p1,p2)
{
  //// find the initial guess for the interior vertices
  // Note: interior vertex 0 is next to vertex 0
  RefVertex p[2];
  p[0] = p1;
  p[1] = p2;
  DVector pv[2];
  DVector grad[2];

  int i;
  for (i=0; i<2; i++) {
      pv[i] = p[i]->point();
      grad[i] = p[i]->gradient();
    }

  for (i=0; i<2; i++) {
      DVector interpv;
      interpv = ((2.0*pv[i])+pv[(i==1)?0:1])*(1.0/3.0);
      RefPoint start(new Point(interpv));
      DVector interpgrad;
      interpgrad = ((2.0*grad[i])+grad[(i==1)?0:1])*(1.0/3.0);
      RefPoint newpoint = vol->solve(start,interpgrad,isovalue);
      vol->gradient(newpoint,interpgrad);
      _interiorvertices[i] = new Vertex(newpoint,interpgrad);
    }

}

double
Edge4::length()
{
  fprintf(stderr,"Edge4::length(): not implemented\n");
  abort();
  return 0.0;
}

/////////////////////////////////////////////////////////////////////////
// Triangle

REF_def(Triangle);
ARRAY_def(RefTriangle);
SET_def(RefTriangle);
ARRAYSET_def(RefTriangle);

Triangle::Triangle(RefEdge v1, RefEdge v2, RefEdge v3)
{
  _edges[0] = v1;
  _edges[1] = v2;
  _edges[2] = v3;

  // edge 0 vertex 0 is (r=0,s=0)
  // edge 0 vertex 1 is (r=0,s=1)
  // edge 1 vertex _vertex is (r=1,s=0)
  // edge 1 vertex (_vertex?0:1) is edge 0 vertex 0

  // swap edges 1 and 2 if necessary
  if (_edges[0]->vertex(0) != _edges[1]->vertex(0)
      && _edges[0]->vertex(0) != _edges[1]->vertex(1)) {
      _edges[1] = v3;
      _edges[2] = v2;
    }

  // find which vertex of edge 1 is not in edge 0
  if (_edges[0]->vertex(0) == _edges[1]->vertex(0)
      || _edges[0]->vertex(1) == _edges[1]->vertex(0)) {
      _edge_vertex = 1;
    }
  else _edge_vertex = 0;

  // make sure that the triangle is really a triangle
  if ( (_edges[0]->vertex(1) != _edges[2]->vertex(0)
        && _edges[0]->vertex(1) != _edges[2]->vertex(1))
       ||
       (_edges[1]->vertex(_edge_vertex) != _edges[2]->vertex(0)
        && _edges[1]->vertex(_edge_vertex) != _edges[2]->vertex(1)) ) {
      fprintf(stderr,"Triangle: given edges that don't form a triangle\n");
      abort();
    }

};

Triangle::~Triangle()
{
}

RefVertex
Triangle::vertex(int i)
{
  if (i==0) {
      return _edges[0]->vertex(0);
    }
  if (i==1) {
      return _edges[0]->vertex(1);
    }
  if (i==2) {
      return _edges[1]->vertex(_edge_vertex);
    }
  return 0;
}

double Triangle::area()
{
  double a = _edges[0]->length();
  double b = _edges[1]->length();
  double c = _edges[2]->length();
  double a2 = a*a;
  double b2 = b*b;
  double c2 = c*c;
  return 0.25 * sqrt( 2.0 * (c2*b2 + c2*a2 + a2*b2)
                      - c2*c2 - b2*b2 - a2*a2);
}

void Triangle::add_vertices(SetRefVertex&set)
{
  set.add(_edges[0]->vertex(0));
  set.add(_edges[0]->vertex(1));
  set.add(_edges[1]->vertex(0));
  set.add(_edges[1]->vertex(1));
  set.add(_edges[2]->vertex(0));
  set.add(_edges[2]->vertex(1));
}

void Triangle::add_edges(SetRefEdge&set)
{
  set.add(_edges[0]);
  set.add(_edges[1]);
}

double
Triangle::interpolate(double r,double s,RefVertex result)
{
  //  ^
  //  |s
  //  3
  //  |\
  //  | \
  //  |  \
  //  |   \
  //  |    \
  //  1-----2-> r
  //
  // use 3 point interpolation
  RefVertex p1(vertex(0));
  RefVertex p2(vertex(2));
  RefVertex p3(vertex(1));

  double L1 = 1 - r - s;
  double L2 = r;
  double L3 = s;

  // Interpolate the gradient
  result->gradient() = L1 * p1->gradient()
    + L2 * p2->gradient()
    + L3 * p3->gradient();

  // Interpolate the point
  DVector x1(p1->point());
  DVector x2(p2->point());
  DVector x3(p3->point());

  DVector newpoint = L1 * x1 + L2 * x2 + L3 * x3;

  for (int i=0; i<result->point()->dimension(); i++) {
      result->point()->operator[](i) = newpoint[i];
    }

  // Find the surface element
  DVector x21 = x2 - x1;
  DVector x31 = x3 - x1;
  return (x21.cross(x31)).norm();
}

/////////////////////////////////////////////////////////////////////////
// Triangle10

REF_def(Triangle10);
ARRAY_def(RefTriangle10);
SET_def(RefTriangle10);
ARRAYSET_def(RefTriangle10);

Triangle10::~Triangle10()
{
}

Triangle10::Triangle10(RefEdge4 v1, RefEdge4 v2, RefEdge4 v3,
                       RefVolume vol,
                       double isovalue):
  Triangle(v1.pointer(),v2.pointer(),v3.pointer())
{
  // set up the 10 vertex array using the 9 vertices from the edges
  // and the midvertices
  //
  // This diagram (Zauhar, Morgon, J Comp Chem, 9(2)171, 1988)
  // gives the number of the vertex and its coordinate in the
  // (r,s) parametric coordinate system.
  //
  //  8 (0,1)
  //  |\
  //  |  \
  //  |    \
  //  |      \
  //  6(0,2/3) 7(1/3,2/3)
  //  |          \
  //  |            \
  //  |              \
  //  4(0,1/3)  9      5(2/3,1/3)
  //  |       (1/3,1/3)  \
  //  |                    \
  //  |                      \
  //  0--------1-------2------3 (1,0)
  //  (0,0)   (1/3,0)  (2/3,0)
  //

  //// put the edges in the same order that Triangle uses

  // swap edges 1 and 2 if necessary
  // (let v2 define s = 0 and v3 define r+s = 1)
  if (v1->vertex(0) != v2->vertex(0)
      && v1->vertex(0) != v2->vertex(1)) {
      RefEdge4 tmp = v2;
      v2 = v3;
      v3 = tmp;
    }

  ////  Find the vertices due to the edges

  // let v1 define r = 0
  _vertices[0] = v1->vertex(0);
  _vertices[4] = v1->interior_vertex(0);
  _vertices[6] = v1->interior_vertex(1);
  _vertices[8] = v1->vertex(1);

  int vertex3;
  int vertex0;
  if (v2->vertex(0) == _vertices[0]) {
      vertex0 = 0;
      vertex3 = 1;
    }
  else {
      vertex0 = 1;
      vertex3 = 0;
    }
  _vertices[3] = v2->vertex(vertex3);
  _vertices[2] = v2->interior_vertex(vertex3);
  _vertices[1] = v2->interior_vertex(vertex0);

  int vertex8;
  if (v3->vertex(0) == _vertices[8]) {
      vertex8 = 0;
      vertex3 = 1;
    }
  else {
      vertex8 = 1;
      vertex3 = 0;
    }
  _vertices[5] = v3->interior_vertex(vertex3);
  _vertices[7] = v3->interior_vertex(vertex8);

  //// find the midvertex of the triangle

  DVector p0(_vertices[0]->point());
  DVector p3(_vertices[3]->point());
  DVector p8(_vertices[8]->point());
  // this vertex corresponds to r = s = 1/3 on the flat triangle
  DVector p9 = (1.0/3.0)*(p0 + p3 + p8);
  RefPoint trialp9(new Point(p9));

  //// find the approximate gradient at the midvertex of the triangle

  DVector& grad0 = _vertices[0]->gradient();
  DVector& grad3 = _vertices[3]->gradient();
  DVector& grad8 = _vertices[8]->gradient();
  DVector grad9 = (1.0/3.0)*(grad0 + grad3 + grad8);

  //// put the midvertex on the surface defined by vol(p9) = isovalue

  RefPoint newpoint = vol->solve(trialp9,grad9,isovalue);
  // compute the true gradient
  vol->gradient(newpoint,grad9);
  _vertices[9] = new Vertex(newpoint,grad9);
};

double Triangle10::area()
{
  fprintf(stderr,"Triangle10::area not implemented\n");
  abort();
}

double
Triangle10::interpolate(double r,double s,RefVertex result)
{
  int i;

  int dim = _vertices[0]->point()->dimension();

  if (result->point()->dimension() != dim) {
      result->point()->resize(dim);
    }
  if (result->gradient().dim() != dim) {
      result->gradient().resize(dim);
    }

  // use 10 point interpolation

  double L1 = 1 - r - s;
  double L2 = r;
  double L3 = s;

  double L1s = -1.0;
  double L2s = 0.0;
  double L3s = 1.0;

  double L1r = -1.0;
  double L2r = 1.0;
  double L3r = 0.0;

  double L1_31 = 3 * L1 - 1;
  double L1_32 = L1_31 - 1;
  double L2_31 = 3 * L2 - 1;
  double L2_32 = L2_31 - 1;
  double L3_31 = 3 * L3 - 1;
  double L3_32 = L3_31 - 1;

  double L1_31s = 3 * L1s;
  double L1_32s = L1_31s;
  double L2_31s = 3 * L2s;
  double L2_32s = L2_31s;
  double L3_31s = 3 * L3s;
  double L3_32s = L3_31s;

  double L1_31r = 3 * L1r;
  double L1_32r = L1_31r;
  double L2_31r = 3 * L2r;
  double L2_32r = L2_31r;
  double L3_31r = 3 * L3r;
  double L3_32r = L3_31r;

  // the interpolation coefficients
  double N[10];
  N[0] = 0.5 * L1 * L1_31 * L1_32;
  N[1] = 4.5 * L1 * L2 * L1_31;
  N[2] = 4.5 * L1 * L2 * L2_31;
  N[3] = 0.5 * L2 * L2_31 * L2_32;
  N[4] = 4.5 * L1 * L3 * L1_31;
  N[5] = 4.5 * L2 * L3 * L2_31;
  N[6] = 4.5 * L1 * L3 * L3_31;
  N[7] = 4.5 * L2 * L3 * L3_31;
  N[8] = 0.5 * L3 * L3_31 * L3_32;
  N[9] = 27  * L1 * L2 * L3;

  // the r derivatives of the interpolation coefficients
  double Nr[10];
  Nr[0] = 0.5*(L1r*L1_31*L1_32 + L1*L1_31r*L1_32 + L1*L1_31*L1_32r);
  Nr[1] = 4.5*(L1r*L2*L1_31 + L1*L2r*L1_31 + L1*L2*L1_31r);
  Nr[2] = 4.5*(L1r*L2*L2_31 + L1*L2r*L2_31 + L1*L2*L2_31r);
  Nr[3] = 0.5*(L2r*L2_31*L2_32 + L2*L2_31r*L2_32 + L2*L2_31*L2_32r);
  Nr[4] = 4.5*(L1r*L3*L1_31 + L1*L3r*L1_31 + L1*L3*L1_31r);
  Nr[5] = 4.5*(L2r*L3*L2_31 + L2*L3r*L2_31 + L2*L3*L2_31r);
  Nr[6] = 4.5*(L1r*L3*L3_31 + L1*L3r*L3_31 + L1*L3*L3_31r);
  Nr[7] = 4.5*(L2r*L3*L3_31 + L2*L3r*L3_31 + L2*L3*L3_31r);
  Nr[8] = 0.5*(L3r*L3_31*L3_32 + L3*L3_31r*L3_32 + L3*L3_31*L3_32r);
  Nr[9] = 27 *(L1r*L2*L3 + L1*L2r*L3 + L1*L2*L3r);

  // the s derivatives of the interpolation coefficients
  double Ns[10];
  Ns[0] = 0.5*(L1s*L1_31*L1_32 + L1*L1_31s*L1_32 + L1*L1_31*L1_32s);
  Ns[1] = 4.5*(L1s*L2*L1_31 + L1*L2s*L1_31 + L1*L2*L1_31s);
  Ns[2] = 4.5*(L1s*L2*L2_31 + L1*L2s*L2_31 + L1*L2*L2_31s);
  Ns[3] = 0.5*(L2s*L2_31*L2_32 + L2*L2_31s*L2_32 + L2*L2_31*L2_32s);
  Ns[4] = 4.5*(L1s*L3*L1_31 + L1*L3s*L1_31 + L1*L3*L1_31s);
  Ns[5] = 4.5*(L2s*L3*L2_31 + L2*L3s*L2_31 + L2*L3*L2_31s);
  Ns[6] = 4.5*(L1s*L3*L3_31 + L1*L3s*L3_31 + L1*L3*L3_31s);
  Ns[7] = 4.5*(L2s*L3*L3_31 + L2*L3s*L3_31 + L2*L3*L3_31s);
  Ns[8] = 0.5*(L3s*L3_31*L3_32 + L3*L3_31s*L3_32 + L3*L3_31*L3_32s);
  Ns[9] = 27 *(L1s*L2*L3 + L1*L2s*L3 + L1*L2*L3s);

  // Interpolate the gradient
  result->gradient().zero();
  for (i=0; i<10; i++)
    result->gradient() += N[i] * _vertices[i]->gradient();

  // Interpolate the point
  DVector x[10];
  for (i=0; i<10; i++)
    x[i] = _vertices[i]->point();

  DVector newpoint(_vertices[0]->point()->dimension());
  newpoint.zero();
  for (i=0; i<10; i++)
    newpoint += N[i] * x[i];

  for (i=0; i<result->point()->dimension(); i++) {
      result->point()->operator[](i) = newpoint[i];
    }

  // Find the surface element
  DVector xr(_vertices[0]->point()->dimension());
  xr.zero();
  for (i=0; i<10; i++)
    xr += Nr[i] * x[i];

  DVector xs(_vertices[0]->point()->dimension());
  xs.zero();
  for (i=0; i<10; i++)
    xs += Ns[i] * x[i];

  // _vertices[9]->point()->print();

  double surface_element = (xr.cross(xs)).norm();
  return surface_element;
}

/////////////////////////////////////////////////////////////////////////
// TriangleIntegrator

TriangleIntegrator::TriangleIntegrator(int order):
  _n(order)
{
  _r = new double [_n];
  _s = new double [_n];
  _w = new double [_n];
}

TriangleIntegrator::~TriangleIntegrator()
{
  delete[] _r;
  delete[] _s;
  delete[] _w;
}

void
TriangleIntegrator::set_w(int i,double w)
{
  _w[i] = w;
}

void
TriangleIntegrator::set_r(int i,double r)
{
  _r[i] = r;
}

void
TriangleIntegrator::set_s(int i,double s)
{
  _s[i] = s;
}

/////////////////////////////////////////////////////////////////////////
// GaussTriangleIntegrator

GaussTriangleIntegrator::GaussTriangleIntegrator(int order):
  TriangleIntegrator(order)
{
  if (order == 1) {
      set_r(0, 1.0/3.0);
      set_s(0, 1.0/3.0);
      set_w(0, 1.0);
    }
  else if (order == 3) {
      set_r(0, 1.0/6.0);
      set_r(1, 2.0/3.0);
      set_r(2, 1.0/6.0);
      set_s(0, 1.0/6.0);
      set_s(1, 1.0/6.0);
      set_s(2, 2.0/3.0);
      set_w(0, 1.0/3.0);
      set_w(1, 1.0/3.0);
      set_w(2, 1.0/3.0);
    }
  else if (order == 4) {
      set_r(0, 1.0/3.0);
      set_r(1, 1.0/5.0);
      set_r(2, 3.0/5.0);
      set_r(3, 1.0/5.0);
      set_s(0, 1.0/3.0);
      set_s(1, 1.0/5.0);
      set_s(2, 1.0/5.0);
      set_s(3, 3.0/5.0);
      set_w(0, -27.0/48.0);
      set_w(1, 25.0/48.0);
      set_w(2, 25.0/48.0);
      set_w(3, 25.0/48.0);
    }
  else if (order == 6) {
      set_r(0, 0.091576213509771);
      set_r(1, 0.816847572980459);
      set_r(2, r(0));
      set_r(3, 0.445948490915965);
      set_r(4, 0.108103018168070);
      set_r(5, r(3));
      set_s(0, r(0));
      set_s(1, r(0));
      set_s(2, r(1));
      set_s(3, r(3));
      set_s(4, r(3));
      set_s(5, r(4));
      set_w(0, 0.109951743655322);
      set_w(1, w(0));
      set_w(2, w(0));
      set_w(3, 0.223381589678011);
      set_w(4, w(3));
      set_w(5, w(3));
    }
  else if (order == 7) {
      set_r(0, 1.0/3.0);
      set_r(1, 0.101286507323456);
      set_r(2, 0.797426985353087);
      set_r(3, r(1));
      set_r(4, 0.470142064105115);
      set_r(5, 0.059715871789770);
      set_r(6, r(4));
      set_s(0, r(0));
      set_s(1, r(1));
      set_s(2, r(1));
      set_s(3, r(2));
      set_s(4, r(4));
      set_s(5, r(4));
      set_s(6, r(5));
      set_w(0, 0.225);
      set_w(1, 0.125939180544827);
      set_w(2, w(1));
      set_w(3, w(1));
      set_w(4, 0.132394152788506);
      set_w(5, w(4));
      set_w(6, w(4));
    }
  else {
      fprintf(stderr,"GaussTriangleIntegrator: invalid order %d\n",order);
      abort();
    }

  // scale the weights by the area of the unit triangle
  for (int i=0; i<order; i++) {
      set_w(i, w(i)*0.5);
    }
}

GaussTriangleIntegrator::~GaussTriangleIntegrator()
{
}

/////////////////////////////////////////////////////////////////////////
// TriangulatedSurface

TriangulatedSurface::TriangulatedSurface():
  _triangle_vertex(0),
  _triangle_edge(0),
  _edge_vertex(0),
  _integrator(new GaussTriangleIntegrator(1))
{
  clear();
}

TriangulatedSurface::~TriangulatedSurface()
{
  clear();
  delete _integrator;
}

void
TriangulatedSurface::set_integrator(TriangleIntegrator*i)
{
  delete _integrator;
  _integrator = i;
}

TriangleIntegrator*
TriangulatedSurface::integrator(int itri)
{
  // currently itri is ignored
  return _integrator;
}

void
TriangulatedSurface::clear()
{
  _completed_surface = 0;

  if (_triangle_vertex) {
      for (int i=0; i<_triangles.length(); i++) {
          delete[] _triangle_vertex[i];
        }
      delete[] _triangle_vertex;
    }
  _triangle_vertex = 0;

  if (_triangle_edge) {
      for (int i=0; i<_triangles.length(); i++) {
          delete[] _triangle_edge[i];
        }
      delete[] _triangle_edge;
    }
  _triangle_edge = 0;

  if (_edge_vertex) {
      for (int i=0; i<_edges.length(); i++) {
          delete[] _edge_vertex[i];
        }
      delete[] _edge_vertex;
    }
  _edge_vertex = 0;

  _have_values = 0;
  _values.clear();

  _vertices.clear();
  _edges.clear();
  _triangles.clear();
}

void
TriangulatedSurface::remove_short_edges(double length_cutoff)
{
  if (!_completed_surface) {
      complete_surface();
    }
  
  int nvertex = _vertices.length();
  int nedge = _edges.length();
  int ntriangle = _triangles.length();

  printf("TriangulatedSurface::remove_short_edges:\n");
  printf("  nvertex nedge ntriangle\n");
  printf("    %3d    %3d     %3d     (%s)\n",
         nvertex,nedge,ntriangle,"initial");

  int* vertex_map = new int[nvertex];
  int* edge_map = new int[nedge];
  int* triangle_map = new int[ntriangle];
  
  for (int ie=0; ie<nedge; ie++) {
      // the edge with potentially two deleted points
      RefEdge E2D(_edges[ie]);
      if (E2D->length() < length_cutoff) {
          int i;
          
          // find the deleted vertices
          ArraysetRefVertex D;
          D.add(_edges[ie]->vertex(0));
          D.add(_edges[ie]->vertex(1));

          // the edges with one deleted point which border triangles with
          // one deleted point.
          ArraysetRefEdge E1D_T1D;
          // the edges with one deleted point which border the triangles
          // with two deleted points.
          ArraysetRefEdge E1D_T2DA;
          ArraysetRefEdge E1D_T2DB;

          // the triangles with two deleted points
          RefTriangle T2DA;
          RefTriangle T2DB;
          // the triangles with one deleted point
          ArraysetRefTriangle T1D;

          // find the triangles T2DA, T2DB, and the triangle set T1D
          // also find the edge arrays E1D_T2DA and E1D_T2DB
          for (i=0; i<ntriangle; i++) {
              if (_triangles[i]->edge(0) == E2D) {
                  if (T2DA.null()) {
                      T2DA = _triangles[i];
                      E1D_T2DA.add(_triangles[i]->edge(1));
                      E1D_T2DA.add(_triangles[i]->edge(2));
                    }
                  else {
                      T2DB = _triangles[i];
                      E1D_T2DB.add(_triangles[i]->edge(1));
                      E1D_T2DB.add(_triangles[i]->edge(2));
                    }
                }
              else if (_triangles[i]->edge(1) == E2D) {
                  if (T2DA.null()) {
                      T2DA = _triangles[i];
                      E1D_T2DA.add(_triangles[i]->edge(0));
                      E1D_T2DA.add(_triangles[i]->edge(2));
                    }
                  else {
                      T2DB = _triangles[i];
                      E1D_T2DB.add(_triangles[i]->edge(0));
                      E1D_T2DB.add(_triangles[i]->edge(2));
                    }
                }
              else if (_triangles[i]->edge(2) == E2D) {
                  if (T2DA.null()) {
                      T2DA = _triangles[i];
                      E1D_T2DA.add(_triangles[i]->edge(0));
                      E1D_T2DA.add(_triangles[i]->edge(1));
                    }
                  else {
                      T2DB = _triangles[i];
                      E1D_T2DB.add(_triangles[i]->edge(0));
                      E1D_T2DB.add(_triangles[i]->edge(1));
                    }
                }
              else if (_triangles[i]->vertex(0) == D[0]
                       || _triangles[i]->vertex(0) == D[1]
                       || _triangles[i]->vertex(1) == D[0]
                       || _triangles[i]->vertex(1) == D[1]
                       || _triangles[i]->vertex(2) == D[0]
                       || _triangles[i]->vertex(2) == D[1]) {
                  T1D.add(_triangles[i]);
                }
            }

          // find the edge array E1D_T1D
          for (i=0; i<T1D.length(); i++) {
              for (int j=0; j<3; j++) {
                  if (T1D[i]->edge(j)->vertex(0) == D[0]
                      ||T1D[i]->edge(j)->vertex(0) == D[1]
                      ||T1D[i]->edge(j)->vertex(1) == D[0]
                      ||T1D[i]->edge(j)->vertex(1) == D[1]) {
                      E1D_T1D.add(T1D[i]->edge(j));
                    }
                }
            }

          // find the vertex replacement, D_rep
          RefVertex D_rep(new Vertex(new Point(3),*(new DVector(3))));
          for (i=0; i<3; i++) {
              D_rep->point()->operator[](i) =
                0.5*(D[0]->point()->operator[](i)
                     + D[1]->point()->operator[](i));
              D_rep->gradient().operator[](i) =
                0.5*(D[0]->gradient().operator[](i)
                     + D[1]->gradient().operator[](i));
            }

          // find the edge replacements, E1D_T2DA_rep
          RefEdge E1D_T2DA_rep;
          if (D[0] != E1D_T2DA[0]->vertex(0)
              && D[1] != E1D_T2DA[0]->vertex(0)) {
              E1D_T2DA_rep = new Edge(E1D_T2DA[0]->vertex(0),D_rep);
            }
          else {
              E1D_T2DA_rep = new Edge(E1D_T2DA[0]->vertex(1),D_rep);
            }

          // find the edge replacements, E1D_T2DB_rep
          RefEdge E1D_T2DB_rep;
          if (D[0] != E1D_T2DB[0]->vertex(0)
              && D[1] != E1D_T2DB[0]->vertex(0)) {
              E1D_T2DB_rep = new Edge(E1D_T2DB[0]->vertex(0),D_rep);
            }
          else {
              E1D_T2DB_rep = new Edge(E1D_T2DB[0]->vertex(1),D_rep);
            }

          // find the edge replacements, E1D_T1D_rep
          ArrayRefEdge E1D_T1D_rep(E1D_T1D.length());
          for (i=0; i<E1D_T1D.length(); i++) {
              if (E1D_T1D[i]->vertex(0) == D[0]
                  || E1D_T1D[i]->vertex(0) == D[1]) {
                  E1D_T1D_rep[i] = new Edge(D_rep,E1D_T1D[i]->vertex(1));
                }
              else {
                  E1D_T1D_rep[i] = new Edge(D_rep,E1D_T1D[i]->vertex(0));
                }
            }

          // find the triangle replacments, T1D_rep
          ArrayRefTriangle T1D_rep(T1D.length());
          for (i=0; i<T1D.length(); i++) {
              ArrayRefEdge edges(3);
              for (int j=0; j<3; j++) {
                  edges[j] = T1D[i]->edge(j);
                  int iedge = E1D_T1D.iseek(edges[j]);
                  if (iedge >= 0) {
                      edges[j] = E1D_T1D_rep[iedge];
                      continue;
                    }
                  iedge = E1D_T2DA.iseek(edges[j]);
                  if (iedge >= 0) {
                      edges[j] = E1D_T2DA_rep;
                      continue;
                    }
                  iedge = E1D_T2DB.iseek(edges[j]);
                  if (iedge >= 0) {
                      edges[j] = E1D_T2DB_rep;
                      continue;
                    }
                }
              T1D_rep[i] = new Triangle(edges[0],edges[1],edges[2]);
            }

          for (i=0; i<ntriangle; i++) triangle_map[i] = i;
          for (i=0; i<nvertex; i++) vertex_map[i] = i;
          for (i=0; i<nedge; i++) edge_map[i] = i;

          // compute the mapping of old vertex indices to new vertex indices
          for (i=0; i<2; i++) {
              int j;
              int jstart = _vertices.iseek(D[i]);
              vertex_map[jstart] = nvertex-2; // delete 2 add 1
              for (j=jstart+1; j<nvertex; j++) {
                  vertex_map[j]--;
                }
            }

          // compute the mapping of old edge indices to new edge indices
          int istart = _edges.iseek(E2D);
          edge_map[istart] = -1;
          for (i=istart+1; i<nedge; i++) edge_map[i]--;
          for (i=0; i<2; i++) {
              int j;
              int jstart = _edges.iseek(E1D_T2DA[i]);
              edge_map[jstart] = nedge - 5; // delete 5 add 1
              for (j=jstart+1; j<nedge; j++) edge_map[j]--;
              jstart = _edges.iseek(E1D_T2DB[i]);
              edge_map[jstart] = nedge - 4; // delete 5 add 2
              for (j=jstart+1; j<nedge; j++) edge_map[j]--;
            }

          // compute the mapping of old triangle indices to new
          istart = _triangles.iseek(T2DA);
          triangle_map[istart] = -1;
          for (i=istart+1; i<ntriangle; i++) triangle_map[i]--;
          istart = _triangles.iseek(T2DB);
          triangle_map[istart] = -1;
          for (i=istart+1; i<ntriangle; i++) triangle_map[i]--;
          
          // replace the E1D_T1D edges
          for (i=0; i<E1D_T1D.length(); i++) {
              _edges[_edges.iseek(E1D_T1D[i])] = E1D_T1D_rep[i];
            }

          // replace the T1D triangles
          for (i=0; i<T1D.length(); i++) {
              _triangles[_triangles.iseek(T1D[i])] = T1D_rep[i];
            }

          // delete the unused objects
          _triangles.del(T2DA);
          _triangles.del(T2DB);
          _edges.del(E2D);
          _vertices.del(D[0]);
          _vertices.del(D[1]);
          for (i=0; i<2; i++) {
              _edges.del(E1D_T2DA[i]);
              _edges.del(E1D_T2DB[i]);
            }

          // put the new vertex on the end
          _vertices.add(D_rep);

          // put the new edges on the end
          _edges.add(E1D_T2DA_rep);
          _edges.add(E1D_T2DB_rep);

          if (_triangle_vertex) {
              int ii;
              for (i=ii=0; ii<ntriangle; ii++) {
                  if (triangle_map[ii] < 0) {
                      delete[] _triangle_vertex[ii];
                      continue;
                    }
                  _triangle_vertex[i] = _triangle_vertex[ii];
                  for (int j=0; j<3; j++) {
                      _triangle_vertex[i][j]
                        = vertex_map[_triangle_vertex[i][j]];
                    }
                  i++;
                }
            }
          if (_triangle_edge) {
              int ii;
              for (i=ii=0; ii<ntriangle; ii++) {
                  if (triangle_map[ii] < 0) {
                      delete[] _triangle_edge[ii];
                      continue;
                    }
                  _triangle_edge[i] = _triangle_edge[ii];
                  for (int j=0; j<3; j++) {
                      _triangle_edge[i][j]
                        = edge_map[_triangle_edge[i][j]];
                    }
                  i++;
                }
            }
          if (_edge_vertex) {
              int ii;
              for (i=ii=0; ii<nedge; ii++) {
                  if (edge_map[ii] < 0) {
                      delete[] _edge_vertex[ii];
                      continue;
                    }
                  _edge_vertex[i] = _edge_vertex[ii];
                  for (int j=0; j<2; j++) {
                      _edge_vertex[i][j]
                        = vertex_map[_edge_vertex[i][j]];
                    }
                  i++;
                }
            }
          ntriangle -= 2;
          nedge -= 3;
          nvertex -= 1;
          
        } /* if below cutoff */
    }
  delete[] edge_map;
  delete[] vertex_map;

  printf("    %3d    %3d     %3d     (%s)\n",
         nvertex,nedge,ntriangle,"final");

  printf("    %3d    %3d     %3d     (%s)\n",
         _vertices.length(),_edges.length(),_triangles.length(),"actual");
}

void
TriangulatedSurface::complete_surface()
{
  int i;
  int ntri = ntriangle();
  for (i=0; i<ntri; i++) {
      RefTriangle tri = triangle(i);
      _edges.add(tri->edge(0));
      _edges.add(tri->edge(1));
      _edges.add(tri->edge(2));
    }
  int ne = nedge();
  // The vertices have possibly already been set up.  If so, then
  // they all better be there.
  if (_vertices.length() == 0) {
      for (i=0; i<ne; i++) {
          RefEdge e = edge(i);
          _vertices.add(e->vertex(0));
          _vertices.add(e->vertex(1));
        }
    }

  if (!_triangle_vertex) {
      // construct the array that converts the triangle number and vertex
      // number within the triangle to the overall vertex number
      _triangle_vertex = new int*[ntri];
      for (i=0; i<ntri; i++) {
          _triangle_vertex[i] = new int[3];
          for (int j=0; j<3; j++) _triangle_vertex[i][j] =
            _vertices.iseek(triangle(i)->vertex(j));
        }
    }

  // construct the array that converts the triangle number and edge number
  // within the triangle to the overall edge number
  _triangle_edge = new int*[ntri];
  for (i=0; i<ntri; i++) {
      _triangle_edge[i] = new int[3];
      for (int j=0; j<3; j++) _triangle_edge[i][j] =
        _edges.iseek(triangle(i)->edge(j));
    }

  // construct the array that converts the edge number and vertex number
  // within the edge to the overall vertex number
  _edge_vertex = new int*[ne];
  for (i=0; i<ne; i++) {
      _edge_vertex[i] = new int[2];
      for (int j=0; j<2; j++)
        _edge_vertex[i][j] = _vertices.iseek(_edges[i]->vertex(j));
    }

  _completed_surface = 1;
}

void
TriangulatedSurface::compute_values(RefVolume&vol)
{
  int n = _vertices.length();
  _values.set_length(n);

  for (int i=0; i<n; i++) {
      vol->SetX(_vertices[i]->point());
      _values[i] = vol->value();
    }
  _have_values = 1;
}

double
TriangulatedSurface::area()
{
  double result = 0.0;
  for (int i=0; i<_triangles.length(); i++) {
      result += _triangles[i]->area();
    }
  return result;
}

double
TriangulatedSurface::volume()
{
  double result = 0.0;
  for (int i=0; i<_triangles.length(); i++) {

      // get the vertices of the triangle
      DVector A(_vertices[triangle_vertex(i,0)]->point());
      DVector B(_vertices[triangle_vertex(i,1)]->point());
      DVector C(_vertices[triangle_vertex(i,2)]->point());

      // project the vertices onto the xy plane
      DVector Axy(A); Axy[2] = 0.0;
      DVector Bxy(B); Bxy[2] = 0.0;
      DVector Cxy(C); Cxy[2] = 0.0;

      // construct the legs of the triangle in the xy plane
      DVector BAxy = Bxy - Axy;
      DVector CAxy = Cxy - Axy;

      // find the lengths of the legs of the triangle in the xy plane
      double baxy = sqrt(BAxy.dot(BAxy));
      double caxy = sqrt(CAxy.dot(CAxy));

      // if one of the legs is of length zero, then there is
      // no contribution from this triangle
      if (baxy < 1.e-16 || caxy < 1.e-16) continue;

      // find the sine of the angle between the legs of the triangle
      // in the xy plane
      double costheta = BAxy.dot(CAxy)/(baxy*caxy);
      double sintheta = sqrt(1.0 - costheta*costheta);

      // the area of the triangle in the xy plane
      double areaxy = 0.5 * baxy * caxy * sintheta;

      // the height of the three corners of the triangle
      // (relative to the z plane)
      double hA = A[2];
      double hB = B[2];
      double hC = C[2];

      // the volume of the space under the triangle
      double volume = areaxy * (hA + (hB + hC - 2.0*hA)/3.0);

      // the sum of the gradients of the triangle's three vertices
      // along the projection axis (z)
      double zgrad
        = _vertices[triangle_vertex(i,0)]->gradient()[2]
        + _vertices[triangle_vertex(i,1)]->gradient()[2]
        + _vertices[triangle_vertex(i,2)]->gradient()[2];

      if (zgrad > 0.0) {
          result += volume;
        }
      else {
          result -= volume;
        }

    }

  // If the volume is negative, then the surface gradients were
  // opposite in sign to the direction assumed.  Flip the sign
  // to fix.
  return fabs(result);
}

void
TriangulatedSurface::add_triangle(RefTriangle&t)
{
  if (_completed_surface) clear();
  else if (_triangle_vertex) {
      fprintf(stderr,"TriangulatedSurface::add_triangle(): "
              "didn't expect vertices\n");
      abort();
    }
  _triangles.add(t);
}

void
TriangulatedSurface::initialize_vertices_triangles(int nvert,int ntri)
{
  if (_completed_surface) clear();
  // construct the array that converts the triangle number and vertex number
  // within the triangle to the overall vertex number
  _triangle_vertex = new int*[ntri];
  for (int i=0; i<ntri; i++) {
      _triangle_vertex[i] = new int[3];
    }
}

void
TriangulatedSurface::add_vertex(RefVertex&v)
{
  _vertices.add(v);
}

void
TriangulatedSurface::add_triangle(RefTriangle&t,
                                  int vertex0,
                                  int vertex1,
                                  int vertex2)
{
  if (!_triangle_vertex) {
      fprintf(stderr,"TriangulatedSurface::add_triangle(): "
              "expected vertices\n");
      abort();
    }
  int i = _triangles.length();
  _triangle_vertex[i][0] = vertex0;
  _triangle_vertex[i][1] = vertex1;
  _triangle_vertex[i][2] = vertex2;
  _triangles.add(t);
}

void
TriangulatedSurface::print(FILE*fp)
{
  fprintf(fp,"TriangulatedSurface:\n");
  int i;

  int np = nvertex();
  fprintf(fp," %3d Vertices:\n",np);
  for (i=0; i<np; i++) {
      RefVertex p = vertex(i);
      int dim = p->dimension();
      fprintf(fp,"  %3d:",i);
      for (int j=0; j<dim; j++) {
          fprintf(fp," % 15.10f", p->point()->operator[](j));
        }
      fprintf(fp,"\n");
    }

  int ne = nedge();
  fprintf(fp," %3d Edges:\n",ne);
  for (i=0; i<ne; i++) {
      RefEdge e = edge(i);
      fprintf(fp,"  %3d: %3d %3d\n",i,
              _vertices.iseek(e->vertex(0)),
              _vertices.iseek(e->vertex(1)));
    }

  int nt = ntriangle();
  fprintf(fp," %3d Triangles:\n",nt);
  for (i=0; i<nt; i++) {
      RefTriangle tri = triangle(i);
      fprintf(fp,"  %3d: %3d %3d %3d\n",i,
              _edges.iseek(tri->edge(0)),
              _edges.iseek(tri->edge(1)),
              _edges.iseek(tri->edge(2)));
    }
}

void
TriangulatedSurface::print_vertices_and_triangles(FILE*fp)
{
  fprintf(fp,"TriangulatedSurface:\n");
  int i;

  int np = nvertex();
  fprintf(fp," %3d Vertices:\n",np);
  for (i=0; i<np; i++) {
      RefVertex p = vertex(i);
      int dim = p->dimension();
      fprintf(fp,"  %3d:",i);
      for (int j=0; j<dim; j++) {
          fprintf(fp," % 15.10f", p->point()->operator[](j));
        }
      fprintf(fp,"\n");
    }

  int nt = ntriangle();
  fprintf(fp," %3d Triangles:\n",nt);
  for (i=0; i<nt; i++) {
      RefTriangle tri = triangle(i);
      fprintf(fp,"  %3d: %3d %3d %3d\n",i,
              _triangle_vertex[i][0],
              _triangle_vertex[i][1],
              _triangle_vertex[i][2]);
    }
}

//////////////////////////////////////////////////////////////////////
// TriangulatedSurface10

TriangulatedSurface10::TriangulatedSurface10(RefVolume&vol,double isovalue):
  _vol(vol),
  _isovalue(isovalue)
{
  // improve the integrator
  set_integrator(new GaussTriangleIntegrator(7));
}

TriangulatedSurface10::~TriangulatedSurface10()
{
}

void
TriangulatedSurface10::complete_surface()
{
  TriangulatedSurface::complete_surface();

  _edges4.set_length(nedge());
  _triangles10.set_length(ntriangle());

  int i;
  for (i=0; i<nedge(); i++) {
      _edges4[i] = new Edge4(vertex(edge_vertex(i,0)),
                             vertex(edge_vertex(i,1)),
                             _vol,_isovalue);
      _edges[i] = _edges4[i].pointer();
    }

  for (i=0; i<ntriangle(); i++) {
      _triangles10[i] = new Triangle10(_edges4[triangle_edge(i,0)],
                                       _edges4[triangle_edge(i,1)],
                                       _edges4[triangle_edge(i,2)],
                                       _vol,_isovalue);
      _triangles[i] = _triangles10[i].pointer();
    }
}

double
TriangulatedSurface10::area()
{
  TriangulatedSurfaceIntegrator tsi(*this);

  double area = 0.0;
  for (tsi = 0; tsi; tsi++) {
      area += tsi.w();
    }
  
  return area;
}

double
TriangulatedSurface10::volume()
{
  TriangulatedSurfaceIntegrator tsi(*this);

  double volume = 0.0;
  RefVertex current = tsi.current();
  for (tsi = 0; tsi; tsi++) {
      DVector norm(current->gradient());
      norm.normalize();
      volume += tsi.w() * norm[0] * current->point()->operator[](0);
    }
  
  return volume;
}

void
TriangulatedSurface10::remove_short_edges(double length_cutoff)
{
  fprintf(stderr,"TriangulatedSurface10::remove_short_edges(): "
          "not implemented\n");
  abort();
}

//////////////////////////////////////////////////////////////////////
// TriangulatedSurfaceIntegrator

TriangulatedSurfaceIntegrator::
  TriangulatedSurfaceIntegrator(TriangulatedSurface&ts)
{
  _ts = &ts;

  _itri = 0;
  _irs = 0;

  if (_ts->nvertex()) {
      int n = _ts->vertex(0)->point()->dimension();
      _current = new Vertex();
      _current->point() = new Point(n);
      _current->gradient().resize(n);
    }
}

TriangulatedSurfaceIntegrator::
  ~TriangulatedSurfaceIntegrator()
{
}

int
TriangulatedSurfaceIntegrator::
  vertex_number(int i)
{
  return _ts->triangle_vertex(_itri,i);
}

RefVertex
TriangulatedSurfaceIntegrator::
  current()
{
  return _current;
}

int
TriangulatedSurfaceIntegrator::
  operator int()
{
  if (_itri < 0 || _itri >= _ts->ntriangle()) return 0;

  TriangleIntegrator* i = _ts->integrator(_itri);
  _s = i->s(_irs);
  _r = i->r(_irs);
  _weight = i->w(_irs);
  _surface_element = _ts->triangle(_itri)->interpolate(_r,_s,_current);

  //printf("%3d: r=%5.3f s=%5.3f w=%8.5f dA=%8.5f ",
  //       _itri, _r, _s, _weight, _surface_element);
  //_current->print(stdout);

  return 1;
}

void
TriangulatedSurfaceIntegrator::
  operator ++()
{
  int n = _ts->integrator(_itri)->n();
  if (_irs == n-1) {
      _irs = 0;
      _itri++;
    }
  else {
      _irs++;
    }
}

int
TriangulatedSurfaceIntegrator::
  operator = (int i)
{
  _itri = i;
  _irs = 0;
  return i;
}


//////////////////////////////////////////////////////////////////////
// UniformLattice

UniformLattice::UniformLattice()
{
  _dim = 0;
  _start = 0;
  _incr = 0;
}


UniformLattice::UniformLattice(int d0,double s0,double i0,
                               int d1,double s1,double i1,
                               int d2,double s2,double i2)
{
  _ndim = 3;
  _dim = new int[3];
  _dim[0] = d0;
  _dim[1] = d1;
  _dim[2] = d2;
  _start = new double[3];
  _incr = new double[3];
  _start[0] = s0;
  _start[1] = s1;
  _start[2] = s2;
  _incr[0] = i0;
  _incr[1] = i1;
  _incr[2] = i2;
}

void
UniformLattice::set_lattice_parameters(int dim,int*dims,
                                       double*start,double*incr)
{
  _ndim = dim;
  _dim = new int[dim];
  _start = new double[dim];
  _incr = new double[dim];
  for (int i=0; i<dim; i++) {
      _dim[i] = dims[i];
      _start[i] = start[i];
      _incr[i] = incr[i];
    }
}

UniformLattice::~UniformLattice()
{
  if (_dim) delete[] _dim;
  if (_start) delete[] _start;
  if (_incr) delete[] _incr;
}

double
UniformLattice::value(int i0,int i1,int i2)
{
  int i[3];
  i[0] = i0;
  i[1] = i1;
  i[2] = i2;
  return value(i);
}

RefPoint
UniformLattice::interpolate(int i0,int i1,int i2,
                            int j0,int j1,int j2,
                            double value)
{
  int i[3];
  i[0] = i0;
  i[1] = i1;
  i[2] = i2;
  int j[3];
  j[0] = j0;
  j[1] = j1;
  j[2] = j2;
  return interpolate(i,j,value);
}

StoredUniformLattice::StoredUniformLattice(int d0,double s0,double i0,
                                           int d1,double s1,double i1,
                                           int d2,double s2,double i2):
  UniformLattice(d0,s0,i0,
                 d1,s1,i1,
                 d2,s2,i2)
{
  _dim = new int[3];
  _dim[0] = d0;
  _dim[1] = d1;
  _dim[2] = d2;
  _data = new double[d0*d1*d2];
  _start = new double[3];
  _incr = new double[3];
  _start[0] = s0;
  _start[1] = s1;
  _start[2] = s2;
  _incr[0] = i0;
  _incr[1] = i1;
  _incr[2] = i2;
}

double
StoredUniformLattice::value(int* i)
{
  return _data[(i[2] + _dim[2]*(i[1] + _dim[1]*(i[0])))];
}

void
StoredUniformLattice::set_value(int i0,int i1,int i2,double v)
{
  _data[(i2 + _dim[2]*(i1 + _dim[1]*(i0)))] = v;
}

StoredUniformLattice::~StoredUniformLattice()
{
  delete[] _data;
}

ImplicitUniformLattice::ImplicitUniformLattice(RefVolume& vol,
                                               int d0,double s0,double i0,
                                               int d1,double s1,double i1,
                                               int d2,double s2,double i2):
  UniformLattice(d0,s0,i0,
                 d1,s1,i1,
                 d2,s2,i2),
  _vol(vol)
{
}

ImplicitUniformLattice::ImplicitUniformLattice(RefVolume& vol,
                                               double resolution,
                                               double minval,
                                               double maxval):
  _vol(vol)
{
  int dim = vol->GetDim();
  int* dims = new int[dim];
  double* start = new double[dim];
  double* incr = new double[dim];
  Point p1(dim),p2(dim);
  vol->boundingbox(minval, maxval, p1, p2);

  int i;
  for (i=0; i<dim; i++) {
      dims[i] = (int)((p2[i]-p1[i])/resolution);
      if (dims[i] < 3) dims[i] = 3;
      incr[i] = (p2[i]-p1[i])/dims[i];
      start[i] = p1[i];
    }

  set_lattice_parameters(dim,dims,start,incr);

  delete[] dims;
  delete[] start;
  delete[] incr;
}

ImplicitUniformLattice::~ImplicitUniformLattice()
{
}

double
ImplicitUniformLattice::value(int* i)
{
  Point p(start(0)+incr(0)*i[0],start(1)+incr(1)*i[1],start(2)+incr(2)*i[2]);
  _vol->SetX(p);
  return _vol->value();
}

RefPoint
ImplicitUniformLattice::interpolate(int* i1,int* i2,double value)
{
  int i;

  RefPoint p1 = new Point(_ndim);
  for (i=0; i<_ndim; i++) p1->operator[](i) = start(i)+incr(i)*i1[i];
  RefPoint p2 = new Point(_ndim);
  for (i=0; i<_ndim; i++) p2->operator[](i) = start(i)+incr(i)*i2[i];

  return _vol->interpolate(p1,p2,value);
}

void
ImplicitUniformLattice::gradient(RefPoint p,DVector& grad)
{
  _vol->gradient(p,grad);
}

IsosurfaceGen::IsosurfaceGen()
{
}

IsosurfaceGen::~IsosurfaceGen()
{
}

MCubesIsosurfaceGen::MCubesIsosurfaceGen(UniformLattice&lat):
  _lattice(lat)
{
}

MCubesIsosurfaceGen::~MCubesIsosurfaceGen()
{
}

void
MCubesIsosurfaceGen::isosurface(double value,
                                TriangulatedSurface& surf)
{
  surf.clear();

  if (_lattice.ndim() != 3) {
      fprintf(stderr,"MCubesIsosurfaceGen::isosurface: "
              " lattice dimension is not 3\n");
      abort();
    }

  int type = P3D_CVTX;
  float* data = new float[_lattice.dim(0)*_lattice.dim(1)*_lattice.dim(2)];
  float* valdata = 0;
  int nx_in = _lattice.dim(0);
  int ny_in = _lattice.dim(1);
  int nz_in = _lattice.dim(2);
  P_Point corner1;
  P_Point corner2;
  int show_inside = 0;
  int ftn_order = 0;

  corner1.x = _lattice.start(0);
  corner1.y = _lattice.start(1);
  corner1.z = _lattice.start(2);
  corner2.x = _lattice.start(0) + _lattice.incr(0)*_lattice.dim(0);
  corner2.y = _lattice.start(1) + _lattice.incr(1)*_lattice.dim(1);
  corner2.z = _lattice.start(2) + _lattice.incr(2)*_lattice.dim(2);

  int ijk = 0;
  for (int i=0; i<nx_in; i++) {
      for (int j=0; j<ny_in; j++) {
          for (int k=0; k<nz_in; k++) {
              data[ijk] = _lattice.value(i,j,k);
              ijk++;
            }
        }
    }

  pg_isosurface_ts(type,data,valdata,
                nx_in,ny_in,nz_in,
                value,&corner1,&corner2,
                show_inside,ftn_order,
                surf,_lattice);
  surf.complete_surface();
}

