
#ifdef __GNUC__
#pragma implementation
#endif

#include <util/keyval/keyval.h>
#include <math/isosurf/triangle.h>
#include <math/scmat/vector3.h>

/////////////////////////////////////////////////////////////////////////
// Triangle

REF_def(Triangle);
ARRAY_def(RefTriangle);
SET_def(RefTriangle);
ARRAYSET_def(RefTriangle);

Triangle::Triangle(const RefEdge& v1, const RefEdge& v2, const RefEdge& v3,
                   unsigned int orientation0):
  _norm(v1->vertex(0)->dimension())
{
  _orientation[0] = orientation0;
  
  _edges[0] = v1;
  _edges[1] = v2;
  _edges[2] = v3;

  // edge 0 vertex _orientation 0 is (r=0,s=0)
  // edge 1 vertex _orientation 1 is (r=0,s=1)
  // edge 2 vertex _orientation 2 is (r=1,s=0)
  // edge 0 vertex (1-_orientation 0) is edge 1 vertex _orientation 1
  // edge 1 vertex (1-_orientation 1) is edge 2 vertex _orientation 2
  // edge 2 vertex (1-_orientation 2) is edge 0 vertex _orientation 0

  RefEdge *e = _edges;
  unsigned int *o = _orientation;

  // swap edges 1 and 2 if necessary
  if (e[0]->vertex(1-o[0]) != e[1]->vertex(0)) {
      if (e[0]->vertex(1-o[0]) != e[1]->vertex(1)) {
          e[1] = v3;
          e[2] = v2;
        }
    }

  // compute the orientation of _edge[1]
  if (e[0]->vertex(1-o[0]) == e[1]->vertex(0)) {
      o[1] = 0;
    }
  else {
      o[1] = 1;
    }

  // compute the orientation of _edge[2]
  if (e[1]->vertex(1-o[1]) == e[2]->vertex(0)) {
      o[2] = 0;
    }
  else {
      o[2] = 1;
    }

  // make sure that the triangle is really a triangle
  if ( e[0]->vertex(1-o[0]) != e[1]->vertex(o[1])
       || e[1]->vertex(1-o[1]) != e[2]->vertex(o[2])
       || e[2]->vertex(1-o[2]) != e[0]->vertex(o[0]))
    {
      fprintf(stderr,"Triangle: given edges that don't form a triangle\n");
      abort();
    }

  // compute the normal to the surface
  SCVector3 BA = vertex(1)->point() - vertex(0)->point();
  SCVector3 CA = vertex(2)->point() - vertex(0)->point();
  SCVector3 N = BA.cross(CA);
  double n = N.norm();
  if (n < 1.0e-15) {
      _norm.assign(0.0);
    }
  else {
      n = 1.0/n;
      for (int i=0; i<3; i++) {
          _norm[i] = N[i]*n;
        }
    }
};

Triangle::~Triangle()
{
}

RefVertex
Triangle::vertex(int i)
{
  return _edges[i]->vertex(_orientation[i]);
}

unsigned int
Triangle::orientation(const RefEdge& e) const
{
  if (e == _edges[0]) return orientation(0);
  if (e == _edges[1]) return orientation(1);
  return orientation(2);
}

int
Triangle::contains(const RefEdge& e) const
{
  if (_edges[0] == e) return 1;
  if (_edges[1] == e) return 1;
  if (_edges[2] == e) return 1;
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
  for (int i=0; i<3; i++) set.add(_edges[i]->vertex(_orientation[i]));
}

void Triangle::add_edges(SetRefEdge&set)
{
  for (int i=0; i<3; i++) set.add(_edges[i]);
}

double
Triangle::interpolate(double r,double s,RefVertex&result)
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
  RefSCVector x1(p1->point());
  RefSCVector x2(p2->point());
  RefSCVector x3(p3->point());

  RefSCVector newpoint = L1 * x1 + L2 * x2 + L3 * x3;

  result->point().assign(newpoint);

  // Find the surface element
  SCVector3 x21 = x2 - x1;
  SCVector3 x31 = x3 - x1;
  return (x21.cross(x31)).norm();
}

void
Triangle::normal(double r,double s,const RefSCVector&v)
{
  v.assign(_norm);
}

void
Triangle::flip()
{
  int i;
  for (i=0; i<3; i++) {
      _orientation[i] = 1 - _orientation[i];
    }
  _norm = -1.0 * _norm;
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

Triangle10::Triangle10(const RefEdge4& v1,
                       const RefEdge4& av2,
                       const RefEdge4& av3,
                       const RefVolume& vol,
                       double isovalue, unsigned int orientation0):
  Triangle(v1.pointer(),av2.pointer(),av3.pointer(),orientation0)
{
  // Refs to edges 2 and 3 might get swapped so copy.
  RefEdge4 v2(av2);
  RefEdge4 v3(av3);
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
  if (v2.pointer() != _edges[1].pointer()) {
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

  RefSCVector p0(_vertices[0]->point());
  RefSCVector p3(_vertices[3]->point());
  RefSCVector p8(_vertices[8]->point());
  // this vertex corresponds to r = s = 1/3 on the flat triangle
  RefSCVector p9 = (1.0/3.0)*(p0 + p3 + p8);
  RefSCVector trialp9(p9.dim());
  trialp9.assign(p9);

  //// find the approximate gradient at the midvertex of the triangle

  RefSCVector grad0 = _vertices[0]->gradient();
  RefSCVector grad3 = _vertices[3]->gradient();
  RefSCVector grad8 = _vertices[8]->gradient();
  RefSCVector grad9 = (1.0/3.0)*(grad0 + grad3 + grad8);

  //// put the midvertex on the surface defined by vol(p9) = isovalue

  RefSCVector newpoint = vol->solve(trialp9,grad9,isovalue);
  // compute the true gradient
  vol->set_x(newpoint);
  grad9 = vol->gradient().copy();
  _vertices[9] = new Vertex(newpoint,grad9);
};

double Triangle10::area()
{
  fprintf(stderr,"Triangle10::area not implemented\n");
  abort();
  return 0.0;
}

double
Triangle10::interpolate(double r,double s,RefVertex&result)
{
  int i;

  RefSCDimension dim = _vertices[0]->point().dim();

  if (!(result->point().dim() == dim) || !(result->gradient().dim() == dim)) {
      fprintf(stderr,"Triangle10::interpolate(): need to allocate vertex\n");
      abort();
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
  result->gradient().assign(0.0);
  for (i=0; i<10; i++)
    result->gradient() = result->gradient() + N[i] * _vertices[i]->gradient();

  // Interpolate the point
  RefSCVector x[10];
  for (i=0; i<10; i++)
    x[i] = _vertices[i]->point();

  RefSCVector newpoint(_vertices[0]->point().dim());
  newpoint.assign(0.0);
  for (i=0; i<10; i++)
    newpoint = newpoint + N[i] * x[i];

  result->point().assign(newpoint);

  // Find the surface element
  RefSCVector xr(_vertices[0]->point().dim());
  xr.assign(0.0);
  for (i=0; i<10; i++)
    xr = xr + Nr[i] * x[i];

  RefSCVector xs(_vertices[0]->point().dim());
  xs.assign(0.0);
  for (i=0; i<10; i++)
    xs = xs + Ns[i] * x[i];

  // _vertices[9]->point()->print();

  SCVector3 xr3(xr);
  SCVector3 xs3(xs);
  double surface_element = (xr3.cross(xs3)).norm();
  return surface_element;
}

/////////////////////////////////////////////////////////////////////////
// TriangleIntegrator

#define CLASSNAME TriangleIntegrator
#define PARENTS public DescribedClass
#define HAVE_KEYVAL_CTOR
//#include <util/state/statei.h>
#include <util/class/classi.h>
void *
TriangleIntegrator::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

TriangleIntegrator::TriangleIntegrator(int order):
  _n(order)
{
  _r = new double [_n];
  _s = new double [_n];
  _w = new double [_n];
}

TriangleIntegrator::TriangleIntegrator(const RefKeyVal& keyval)
{
  _n = keyval->intvalue("n");
  if (keyval->error() != KeyVal::OK) {
      _n = 7;
    }
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
TriangleIntegrator::set_n(int n)
{
  delete[] _r;
  delete[] _s;
  delete[] _w;
  _n = n;
  _r = new double [_n];
  _s = new double [_n];
  _w = new double [_n];
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

#define CLASSNAME GaussTriangleIntegrator
#define PARENTS public TriangleIntegrator
#define HAVE_KEYVAL_CTOR
//#include <util/state/statei.h>
#include <util/class/classi.h>
void *
GaussTriangleIntegrator::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = TriangleIntegrator::_castdown(cd);
  return do_castdowns(casts,cd);
}

GaussTriangleIntegrator::GaussTriangleIntegrator(const RefKeyVal& keyval):
  TriangleIntegrator(keyval)
{
  init_rw(n());
}

GaussTriangleIntegrator::GaussTriangleIntegrator(int order):
  TriangleIntegrator(order)
{
  init_rw(n());
}

void
GaussTriangleIntegrator::set_n(int n)
{
  TriangleIntegrator::set_n(n);
  init_rw(n);
}

void
GaussTriangleIntegrator::init_rw(int order)
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
