
#ifdef __GNUC__
#pragma implementation
#endif

#include <util/keyval/keyval.h>
#include <math/isosurf/triangle.h>
#include <math/scmat/vector3.h>

static inline int
ijk_to_index(int i, int j, int k)
{
  int n = i + j + k;
  int ir = n - i;
  return (ir*(ir+1)>>1) + j;
}

static inline int
order_to_nvertex(int order)
{
  return ((order+1)*(order+2)>>1);
}

/////////////////////////////////////////////////////////////////////////
// Triangle

// Here is the layout of the vertices and edges in the triangles:
//
//                       vertex(1) (r=0, s=1)
//                          +
//                         / \  _edges[1].vertex(_orientation1)
//                        /   \
//                       /     \
//                      /       \
//                     /         \
//                    /           \
//         _edges[0] /             \  _edges[1]
//          (r = 0) /               \   (1-r-s = 0)
//                 /                 \
//                /                   \
//               /                     \
//              /                       \ _edges[1].vertex(!_orientation1)
//             /                         \
//   vertex(0)+---------------------------+
// (r=0, s=0)        _edges[2] (s = 0)      vertex(2) (r=1, s=0)
//
//  Zienkiewicz and Taylor, "The Finite Element Method", 4th Ed, Vol 1,
//  use
//      L1 = 1 - r - s
//      L2 = r,
//      L3 = s.
//  I also use these below.
//

REF_def(Triangle);
ARRAY_def(RefTriangle);
SET_def(RefTriangle);
ARRAYSET_def(RefTriangle);

Triangle::Triangle(const RefEdge& v1, const RefEdge& v2, const RefEdge& v3,
                   unsigned int orientation0)
{
  _orientation0 = orientation0;

  _edges[0] = v1;
  _edges[1] = v2;
  _edges[2] = v3;

  // edge 0 corresponds to r = 0
  // edge 1 corresponds to (1-r-s) = 0
  // edge 2 corresponds to s = 0
  // edge 0 vertex _orientation0 is (r=0,s=0)
  // edge 1 vertex _orientation1 is (r=0,s=1)
  // edge 2 vertex _orientation2 is (r=1,s=0)
  // edge 0 vertex (1-_orientation0) is edge 1 vertex _orientation1
  // edge 1 vertex (1-_orientation1) is edge 2 vertex _orientation2
  // edge 2 vertex (1-_orientation2) is edge 0 vertex _orientation0

  RefEdge *e = _edges;

  // swap edges 1 and 2 if necessary
  if (e[0]->vertex(1-_orientation0) != e[1]->vertex(0)) {
      if (e[0]->vertex(1-_orientation0) != e[1]->vertex(1)) {
          e[1] = v3;
          e[2] = v2;
        }
    }

  // compute the orientation of _edge[1]
  if (e[0]->vertex(1-_orientation0) == e[1]->vertex(0)) {
      _orientation1 = 0;
    }
  else {
      _orientation1 = 1;
    }

  // compute the orientation of _edge[2]
  if (e[1]->vertex(1-_orientation1) == e[2]->vertex(0)) {
      _orientation2 = 0;
    }
  else {
      _orientation2 = 1;
    }

  // make sure that the triangle is really a triangle
  if ( e[0]->vertex(1-_orientation0) != e[1]->vertex(_orientation1)
       || e[1]->vertex(1-_orientation1) != e[2]->vertex(_orientation2)
       || e[2]->vertex(1-_orientation2) != e[0]->vertex(_orientation0))
    {
      fprintf(stderr,"Triangle: given edges that don't form a triangle\n");
      abort();
    }

  _order = 1;
  _vertices = new RefVertex[3];

  _vertices[ijk_to_index(_order, 0, 0)] = vertex(0);
  _vertices[ijk_to_index(0, _order, 0)] = vertex(1);
  _vertices[ijk_to_index(0, 0, _order)] = vertex(2);
}

Triangle::~Triangle()
{
  if (_vertices) delete[] _vertices;
}

RefVertex
Triangle::vertex(int i)
{
  return _edges[i]->vertex(orientation(i));
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

double Triangle::flat_area()
{
  double a = _edges[0]->straight_length();
  double b = _edges[1]->straight_length();
  double c = _edges[2]->straight_length();
  double a2 = a*a;
  double b2 = b*b;
  double c2 = c*c;
  return 0.25 * sqrt( 2.0 * (c2*b2 + c2*a2 + a2*b2)
                      - c2*c2 - b2*b2 - a2*a2);
}

void Triangle::add_vertices(SetRefVertex&set)
{
  for (int i=0; i<3; i++) set.add(_edges[i]->vertex(orientation(i)));
}

void Triangle::add_edges(SetRefEdge&set)
{
  for (int i=0; i<3; i++) set.add(_edges[i]);
}

static inline void
init_coef_deriv(double L, int order, double *Lcoef, double *Lcoefderiv)
{
  int i;
  Lcoef[0] = 1.0;
  Lcoefderiv[0] = 0.0;
  double spacing = 1.0/order;
  for (i=1; i<=order; i++) {
      Lcoef[i] = Lcoef[i-1] * (L - (i-1)*spacing)/(i*spacing);
      Lcoefderiv[i] = Lcoefderiv[i-1] * (L - (i-1)*spacing)/(i*spacing)
                      + Lcoef[i-1]/(i*spacing);
    }
}

double
Triangle::interpolate(double r,double s,const RefVertex&result)
{
  int i, j, k;

  double L1 = 1 - r - s;
  double L2 = r;
  double L3 = s;

  double L1coef[max_order+1];
  double L2coef[max_order+1];
  double L3coef[max_order+1];

  double L1coefderiv[max_order+1];
  double L2coefderiv[max_order+1];
  double L3coefderiv[max_order+1];

  // the r derivatives
  double L1coef_r[max_order+1];
  double L2coef_r[max_order+1];
  double L3coef_r[max_order+1];

  // the s derivatives
  double L1coef_s[max_order+1];
  double L2coef_s[max_order+1];
  double L3coef_s[max_order+1];

  init_coef_deriv(L1, _order, L1coef, L1coefderiv);
  init_coef_deriv(L2, _order, L2coef, L2coefderiv);
  init_coef_deriv(L3, _order, L3coef, L3coefderiv);

  // convert into r and s derivatives
  for (i=0; i<=_order; i++) {
      L1coef_r[i] = -L1coefderiv[i];
      L1coef_s[i] = -L1coefderiv[i];
      L2coef_r[i] =  L2coefderiv[i];
      L2coef_s[i] =  0.0;
      L3coef_r[i] =  0.0;
      L3coef_s[i] =  L3coefderiv[i];
    }

  RefSCVector tmp(_vertices[0]->point()->dim());
  RefSCVector x_s(_vertices[0]->point()->dim());
  RefSCVector x_r(_vertices[0]->point()->dim());
  tmp.assign(0.0);
  x_s.assign(0.0);
  x_r.assign(0.0);
  for (i=0; i<=_order; i++) {
      for (j=0; j <= _order-i; j++) {
          k = _order - i - j;
          tmp.accumulate((L1coef[i]*L2coef[j]*L3coef[k])
                         *_vertices[ijk_to_index(i,j,k)]->point());
          x_s.accumulate((L1coef_s[i]*L2coef[j]*L3coef[k]
                          +L1coef[i]*L2coef_s[j]*L3coef[k]
                          +L1coef[i]*L2coef[j]*L3coef_s[k])
                         *_vertices[ijk_to_index(i,j,k)]->point());
          x_r.accumulate((L1coef_r[i]*L2coef[j]*L3coef[k]
                          +L1coef[i]*L2coef_r[j]*L3coef[k]
                          +L1coef[i]*L2coef[j]*L3coef_r[k])
                         *_vertices[ijk_to_index(i,j,k)]->point());
        }
    }
  result->point().assign(tmp);

  if (result->normal().nonnull()) {
      tmp.assign(0.0);
      for (i=0; i<_order; i++) {
          for (j=0; j <= _order-i; j++) {
              k = _order - i - j;
              tmp.accumulate((L1coef[i]*L2coef[j]*L3coef[k])
                             *_vertices[ijk_to_index(i,j,k)]->point());
            }
        }
      result->normal().assign(tmp);
    }

#undef PRINT_COEFFICIENTS
#ifdef PRINT_COEFFICIENTS
  static int print = 0;
  if (!print) {
      print = 1;
      printf(" r = %10.7f s = %10.7f\n", r, s);
      printf(" L1 = %10.7f L2 = %10.7f L3 = %10.7f\n", L1, L2, L3);
      printf(" The interpolation coefficients:\n");
      for (i=0; i<=_order; i++) {
          printf("  %d   %10.7f %10.7f %10.7f\n",
                 i, L1coef[i], L2coef[i], L3coef[i]);
        }
      printf(" The interpolation coefficients' r derivs:\n");
      for (i=0; i<=_order; i++) {
          printf("  %d   %10.7f %10.7f %10.7f\n",
                 i, L1coef_r[i], L2coef_r[i], L3coef_r[i]);
        }
      printf(" The interpolation coefficients' s derivs:\n");
      for (i=0; i<=_order; i++) {
          printf("  %d   %10.7f %10.7f %10.7f\n",
                 i, L1coef_s[i], L2coef_s[i], L3coef_s[i]);
        }

      printf(" The interpolation coefficient products including derivs:\n");
      for (i=0; i<=_order; i++) {
          for (j=0; j<=_order-i; j++) {
              k = _order - i - j;
              int index = ijk_to_index(i,j,k);
              printf("  %2d (%d %d %d)  %10.7f %10.7f %10.7f\n",
                     index, i, j, k,
                     L1coef[i]*L2coef[j]*L3coef[k],
                     (L1coef_r[i]*L2coef[j]*L3coef[k]
                      +L1coef[i]*L2coef_r[j]*L3coef[k]
                      +L1coef[i]*L2coef[j]*L3coef_r[k]),
                     (L1coef_s[i]*L2coef[j]*L3coef[k]
                      +L1coef[i]*L2coef_s[j]*L3coef[k]
                      +L1coef[i]*L2coef[j]*L3coef_s[k])
                     );
            }
        }

      printf(" The corner vertices:\n");
      for (i=0; i<3; i++) {
          printf(" %d: ", i);
          SCVector3 v(vertex(i)->point());
          v.print();
        }

      printf(" The interpolated vertex: ");
      SCVector3 interp(result->point());
      interp.print();
    }
#endif // PRINT_COEFFICIENTS

  // Find the surface element
  SCVector3 xr3(x_r);
  SCVector3 xs3(x_s);
  double surface_element = (xr3.cross(xs3)).norm();
  return surface_element;
}

void
Triangle::flip()
{
  _orientation0 = _orientation0?0:1;
  _orientation1 = _orientation1?0:1;
  _orientation2 = _orientation2?0:1;
}

void
Triangle::set_order(int order, const RefVolume&vol, double isovalue)
{
  if (order > max_order) {
      fprintf(stderr,"Triangle::set_order: max_order = %d but order = %d\n",
              max_order, order);
      abort();
    }

  int i, j, k;

  static int print = 0;
  if (print == 0) {
      print = 1;
      printf("order = %d new order = %d\n", _order, order);
      for (i=0; i<order_to_nvertex(_order); i++) {
          printf("  _vertices[%d] = 0x%08x\n", i, _vertices[i].pointer());
        }
      for (i=0; i<3; i++) {
          printf("  vertex(%d) = 0x%08x\n", i, vertex(i).pointer());
        }
    }

  if (edge(0)->order() != order
      ||edge(1)->order() != order
      ||edge(2)->order() != order) {
      fprintf(stderr,"Triangle::set_order: edge order doesn't match\n");
      abort();
    }

  _order = order;
  delete[] _vertices;

  _vertices = new RefVertex[order_to_nvertex(_order)];

  // fill in the corner vertices
  _vertices[ijk_to_index(_order, 0, 0)] = vertex(0);
  _vertices[ijk_to_index(0, 0, _order)] = vertex(1);
  _vertices[ijk_to_index(0, _order, 0)] = vertex(2);

  // fill in the interior edge vertices
  for (i = 1; i < _order; i++) {
      j = _order - i;
      _vertices[ijk_to_index(0, i, j)] = _edges[1]->interior_vertex(
          _orientation1?i:j);
      _vertices[ijk_to_index(j, 0, i)] = _edges[0]->interior_vertex(
          _orientation0?i:j);
      _vertices[ijk_to_index(i, j, 0)] = _edges[2]->interior_vertex(
          _orientation2?i:j);
    }

  if (order != _order) abort();

  // find the triangle's interior vertices
  RefSCVector p0(vertex(0)->point());
  RefSCVector p1(vertex(1)->point());
  RefSCVector p2(vertex(2)->point());
  RefSCVector norm0(vertex(0)->normal());
  RefSCVector norm1(vertex(1)->normal());
  RefSCVector norm2(vertex(2)->normal());

  if (order != _order) abort();

  for (i=0; i<=_order; i++) {
      double I = (1.0*i)/_order;
      for (j=0; j<=_order-i; j++) {
          RefSCVector trialpoint;
          RefSCVector trialnorm;
          RefSCVector newpoint;
          double J = (1.0*j)/_order;
          k = _order - i - j;
          if (!i || !j || !k) continue; // interior point check
          double K = (1.0*k)/_order;
          int index = ijk_to_index(i,j,k);
          // this get approximate vertices and normals
          trialpoint = I*p0 + J*p1 + K*p2;
          trialnorm = I*norm0 + J*norm1 + K*norm2;
          // now refine that guess
          newpoint = vol->solve(trialpoint,trialnorm,isovalue);
          // compute the true normal
          vol->set_x(newpoint);
          if (vol->gradient_implemented()) {
              trialnorm = vol->gradient().copy();
            }
          trialnorm.normalize();
          _vertices[index] = new Vertex(newpoint,trialnorm);
        }
    }

  if (print == 1) {
      print = 2;
      for (i=0; i<order_to_nvertex(_order); i++) {
          printf("  _vertices[%d] = 0x%08x\n", i, _vertices[i].pointer());
        }
      for (i=0; i<3; i++) {
          printf("  vertex(%d) = 0x%08x\n", i, vertex(i).pointer());
        }
    }
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
  printf("Created a GaussTriangleIntegrator with n = %d\n", n());
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
