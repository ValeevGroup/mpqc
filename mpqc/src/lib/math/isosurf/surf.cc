
#ifdef __GNUC__
#pragma implementation
#endif

#include <util/keyval/keyval.h>
#include <math/scmat/matrix.h>
#include <math/scmat/vector3.h>
#include <math/isosurf/surf.h>
#include <math/isosurf/isosurf.h>

/////////////////////////////////////////////////////////////////////////
// TriangulatedSurface

#define CLASSNAME TriangulatedSurface
#define PARENTS public DescribedClass
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
//#include <util/state/statei.h>
#include <util/class/classi.h>
void *
TriangulatedSurface::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

TriangulatedSurface::TriangulatedSurface():
  _triangle_vertex(0),
  _triangle_edge(0),
  _edge_vertex(0),
  _integrator(new GaussTriangleIntegrator(1)),
  _index_to_vertex(0),
  _index_to_edge(0),
  _index_to_triangle(0),
  _vertex_to_index(0),
  _edge_to_index(0),
  _triangle_to_index(0),
  _tmp_edges(RefEdgeAVLSet())
{
  clear();
}

TriangulatedSurface::TriangulatedSurface(const RefKeyVal&):
  _triangle_vertex(0),
  _triangle_edge(0),
  _edge_vertex(0),
  _integrator(new GaussTriangleIntegrator(1)),
  _index_to_vertex(0),
  _index_to_edge(0),
  _index_to_triangle(0),
  _vertex_to_index(0),
  _edge_to_index(0),
  _triangle_to_index(0),
  _tmp_edges(RefEdgeAVLSet())
{
  clear();
}

TriangulatedSurface::~TriangulatedSurface()
{
  clear();
}

void
TriangulatedSurface::topology_info(FILE*fp)
{
  topology_info(nvertex(), nedge(), ntriangle(), fp);
}

void
TriangulatedSurface::topology_info(int v, int e, int t, FILE* fp)
{
  // Given v vertices i expect 2*v - 4*n_surface triangles
  // and 3*v - 6*n_surface edges
  fprintf(fp, "n_vertex = %d, n_edge = %d, n_triangle = %d:\n",
          v, e, t);
  int nsurf_e = ((3*v - e)%6 == 0)? (3*v - e)/6 : -1;
  int nsurf_t = ((2*v - t)%4 == 0)? (2*v - t)/4 : -1;
  if ((nsurf_e!=-1) && (nsurf_e == nsurf_t)) {
      fprintf(fp,"  this is consistent with n_closed_surface - n_hole = %d\n",
              nsurf_e);
    }
  else {
      fprintf(fp,"  this implies that some surfaces are not closed\n");
    }
}

void
TriangulatedSurface::set_integrator(const RefTriangleIntegrator& i)
{
  _integrator = i;
}

RefTriangleIntegrator
TriangulatedSurface::integrator(int)
{
  // currently the argument, the integer index of the triangle, is ignored
  return _integrator;
}

void
TriangulatedSurface::clear_int_arrays()
{
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

  _completed_surface = 0;
}

void
TriangulatedSurface::clear()
{
  _completed_surface = 0;

  clear_int_arrays();

  _have_values = 0;
  _values.clear();

  _vertices.clear();
  _edges.clear();
  _triangles.clear();

  _tmp_edges.clear();
}

void
TriangulatedSurface::complete_surface()
{
  complete_ref_arrays();
  complete_int_arrays();

  _completed_surface = 1;
}

void
TriangulatedSurface::complete_ref_arrays()
{
  _tmp_edges.clear();

  int i;
  int ntri = ntriangle();
  _edges.clear();
  for (i=0; i<ntri; i++) {
      RefTriangle tri = triangle(i);
      add_edge(tri->edge(0));
      add_edge(tri->edge(1));
      add_edge(tri->edge(2));
    }
  int ne = nedge();
  _vertices.clear();
  for (i=0; i<ne; i++) {
      RefEdge e = edge(i);
      add_vertex(e->vertex(0));
      add_vertex(e->vertex(1));
    }
}

void
TriangulatedSurface::complete_int_arrays()
{
  clear_int_arrays();
  
  int i;
  int ntri = ntriangle();
  int ne = nedge();

  // construct the array that converts the triangle number and vertex
  // number within the triangle to the overall vertex number
  _triangle_vertex = new int*[ntri];
  for (i=0; i<ntri; i++) {
      _triangle_vertex[i] = new int[3];
      for (int j=0; j<3; j++)
          _triangle_vertex[i][j] =
                   _vertex_to_index[_vertices.seek(triangle(i)->vertex(j))];
    }

  // construct the array that converts the triangle number and edge number
  // within the triangle to the overall edge number
  _triangle_edge = new int*[ntri];
  for (i=0; i<ntri; i++) {
      _triangle_edge[i] = new int[3];
      for (int j=0; j<3; j++)
          _triangle_edge[i][j] =
             _edge_to_index[_edges.seek(triangle(i)->edge(j))];
    }

  // construct the array that converts the edge number and vertex number
  // within the edge to the overall vertex number
  _edge_vertex = new int*[ne];
  for (i=0; i<ne; i++) {
      _edge_vertex[i] = new int[2];
      for (int j=0; j<2; j++)
        _edge_vertex[i][j]
            = _vertex_to_index[_vertices.seek(edge(i)->vertex(j))];
    }
}

void
TriangulatedSurface::compute_values(RefVolume&vol)
{
  int n = _vertices.length();
  _values.set_length(n);

  for (int i=0; i<n; i++) {
      vol->set_x(vertex(i)->point());
      _values[i] = vol->value();
    }
  _have_values = 1;
}

double
TriangulatedSurface::area()
{
  double result = 0.0;
  for (Pix i=_triangles.first(); i; _triangles.next(i)) {
      result += _triangles(i)->area();
    }
  return result;
}

double
TriangulatedSurface::volume()
{
  double result = 0.0;
  for (int i=0; i<_triangles.length(); i++) {

      // get the vertices of the triangle
      SCVector3 A(vertex(triangle_vertex(i,0))->point());
      SCVector3 B(vertex(triangle_vertex(i,1))->point());
      SCVector3 C(vertex(triangle_vertex(i,2))->point());

      // project the vertices onto the xy plane
      SCVector3 Axy(A); Axy[2] = 0.0;
      SCVector3 Bxy(B); Bxy[2] = 0.0;
      SCVector3 Cxy(C); Cxy[2] = 0.0;

      // construct the legs of the triangle in the xy plane
      SCVector3 BAxy = Bxy - Axy;
      SCVector3 CAxy = Cxy - Axy;

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

      // the orientation of the triangle along the projection axis (z)
      SCVector3 BA(B-A);
      SCVector3 CA(C-A);
      double z_orientation = BA.cross(CA)[2];

      if (z_orientation > 0.0) {
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
TriangulatedSurface::add_vertex(const RefVertex&t)
{
  int i = _vertices.length();
  RefVertex tnotconst(t);
  Pix ix = _vertices.add(tnotconst);
  if (i != _vertices.length()) {
      _index_to_vertex[i] = ix;
      _vertex_to_index[ix] = i;
    }
  if (_index_to_vertex.length() != _vertex_to_index.length()) {
      fprintf(stderr,"TriangulatedSurface::add_vertex: length mismatch\n");
      abort();
    }
}

void
TriangulatedSurface::add_edge(const RefEdge&t)
{
  int i = _edges.length();
  RefEdge tnotconst(t);
  Pix ix = _edges.add(tnotconst);
  if (i != _edges.length()) {
      _index_to_edge[i] = ix;
      _edge_to_index[ix] = i;
    }
  if (_index_to_edge.length() != _edge_to_index.length()) {
      fprintf(stderr,"TriangulatedSurface::add_edge: length mismatch\n");
      abort();
    }
}

void
TriangulatedSurface::add_triangle(const RefTriangle&t)
{
  if (_completed_surface) clear();
  int i = _triangles.length();
  RefTriangle tnotconst(t);
  Pix ix = _triangles.add(tnotconst);
  if (i != _triangles.length()) {
      _index_to_triangle[i] = ix;
      _triangle_to_index[ix] = i;
    }
}

void
TriangulatedSurface::add_triangle(const RefVertex& v1,
                                  const RefVertex& v2,
                                  const RefVertex& v3)
{
  // Find this triangle's edges if they have already be created
  // for some other triangle.
  RefEdge e0, e1, e2;

  RefVertex v1_not_const(v1);
  RefEdgeAVLSet v1edges;
  v1edges |= _tmp_edges[v1_not_const];

  RefVertex v2_not_const(v2);
  RefEdgeAVLSet v2edges;
  v2edges |= _tmp_edges[v2_not_const];

  Pix ix;
  for (ix = v1edges.first(); ix; v1edges.next(ix)) {
      RefEdge& e = v1edges(ix);
      if (e->vertex(0) == v2 || e->vertex(1) == v2) {
          e0 = e;
        }
      else if (e->vertex(0) == v3 || e->vertex(1) == v3) {
          e2 = e;
        }
    }
  for (ix = v2edges.first(); ix; v1edges.next(ix)) {
      RefEdge& e = v2edges(ix);
      if (e->vertex(0) == v3 || e->vertex(1) == v3) {
          e1 = e;
        }
    }

  RefVertex v3_not_const(v3);
  if (e0.null()) {
      e0 = newEdge(v1,v2);
      _tmp_edges[v1_not_const].add(e0);
      _tmp_edges[v2_not_const].add(e0);
    }
  if (e1.null()) {
      e1 = newEdge(v2,v3);
      _tmp_edges[v2_not_const].add(e1);
      _tmp_edges[v3_not_const].add(e1);
    }
  if (e2.null()) {
      e2 = newEdge(v3,v1);
      _tmp_edges[v3_not_const].add(e2);
      _tmp_edges[v1_not_const].add(e2);
    }
  
  int orientation;
  if (e0->vertex(0) == v1) {
      orientation = 0;
    }
  else {
      orientation = 1;
    }
  
  add_triangle(newTriangle(e0,e1,e2,orientation));
}

// If a user isn't keeping track of edges while add_triangle is being
// used to build the surface, then this can be called to see if an edge
// already exists (at a great performance cost).
RefEdge
TriangulatedSurface::find_edge(const RefVertex& v1, const RefVertex& v2)
{
  Pix i;

  for (i=_triangles.first(); i; _triangles.next(i)) {
      RefTriangle t = _triangles(i);
      RefEdge e1 = t->edge(0);
      RefEdge e2 = t->edge(1);
      RefEdge e3 = t->edge(2);
      if (e1->vertex(0) == v1 && e1->vertex(1) == v2) return e1;
      if (e1->vertex(1) == v1 && e1->vertex(0) == v2) return e1;
      if (e2->vertex(0) == v1 && e2->vertex(1) == v2) return e2;
      if (e2->vertex(1) == v1 && e2->vertex(0) == v2) return e2;
      if (e3->vertex(0) == v1 && e3->vertex(1) == v2) return e3;
      if (e3->vertex(1) == v1 && e3->vertex(0) == v2) return e3;
    }

  return 0;
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
          fprintf(fp," % 15.10f", (double)p->point()[j]);
        }
      fprintf(fp,"\n");
    }

  int ne = nedge();
  fprintf(fp," %3d Edges:\n",ne);
  for (i=0; i<ne; i++) {
      RefEdge e = edge(i);
      fprintf(fp,"  %3d: %3d %3d\n",i,
              _vertex_to_index[_vertices.seek(e->vertex(0))],
              _vertex_to_index[_vertices.seek(e->vertex(1))]);
    }

  int nt = ntriangle();
  fprintf(fp," %3d Triangles:\n",nt);
  for (i=0; i<nt; i++) {
      RefTriangle tri = triangle(i);
      fprintf(fp,"  %3d: %3d %3d %3d\n",i,
              _edge_to_index[_edges.seek(tri->edge(0))],
              _edge_to_index[_edges.seek(tri->edge(1))],
              _edge_to_index[_edges.seek(tri->edge(2))]);
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
          fprintf(fp," % 15.10f", (double)p->point()[j]);
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

void
TriangulatedSurface::print_geomview_format(FILE*fp)
{
  fprintf(fp,"OFF\n");

  fprintf(fp,"%d %d %d\n",nvertex(),ntriangle(),nedge());
  int i;

  int np = nvertex();
  for (i=0; i<np; i++) {
      RefVertex p = vertex(i);
      int dim = p->dimension();
      for (int j=0; j<dim; j++) {
          fprintf(fp," % 15.10f", (double)p->point()[j]);
        }
      fprintf(fp,"\n");
    }

  int nt = ntriangle();
  for (i=0; i<nt; i++) {
      RefTriangle tri = triangle(i);
      fprintf(fp," 3 %3d %3d %3d\n",
              _triangle_vertex[i][0],
              _triangle_vertex[i][1],
              _triangle_vertex[i][2]);
    }
}

void
TriangulatedSurface::recompute_index_maps()
{
  int i;
  Pix I;

  // fix the index maps
  _vertex_to_index.clear();
  _edge_to_index.clear();
  _triangle_to_index.clear();

  _index_to_vertex.clear();
  _index_to_edge.clear();
  _index_to_triangle.clear();

  int ne = _edges.length();
  int nv = _vertices.length();
  int nt = _triangles.length();

  for (i=0, I = _vertices.first(); I; i++, _vertices.next(I)) {
      _vertex_to_index[I] = i;
      _index_to_vertex[i] = I;
    }

  for (i=0, I = _edges.first(); I; i++, _edges.next(I)) {
      _edge_to_index[I] = i;
      _index_to_edge[i] = I;
    }

  for (i=0, I = _triangles.first(); I; i++, _triangles.next(I)) {
      _triangle_to_index[I] = i;
      _index_to_triangle[i] = I;
    }

}

Edge*
TriangulatedSurface::newEdge(const RefVertex& v0, const RefVertex& v1) const
{
  return new Edge(v0,v1);
}

Triangle*
TriangulatedSurface::newTriangle(const RefEdge& e0,
                                 const RefEdge& e1,
                                 const RefEdge& e2,
                                 int orientation) const
{
  return new Triangle(e0,e1,e2,orientation);
}

//////////////////////////////////////////////////////////////////////
// TriangulatedSurface10

#define CLASSNAME TriangulatedSurface10
#define PARENTS public TriangulatedSurface
#define HAVE_KEYVAL_CTOR
//#include <util/state/statei.h>
#include <util/class/classi.h>
void *
TriangulatedSurface10::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = TriangulatedSurface::_castdown(cd);
  return do_castdowns(casts,cd);
}

TriangulatedSurface10::TriangulatedSurface10(const RefVolume&vol,
                                             double isovalue):
  _vol(vol),
  _isovalue(isovalue)
{
  // improve the integrator
  set_integrator(new GaussTriangleIntegrator(7));
}

TriangulatedSurface10::TriangulatedSurface10(const RefKeyVal& keyval)
{
  _vol = keyval->describedclassvalue("volume");
  _isovalue = keyval->doublevalue("value");

  // If an integrator was not specified in the input, then choose
  // a better one than the TriangulatedSurface default.
  if (!keyval->exists("integrator")) {
      set_integrator(new GaussTriangleIntegrator(7));
    }
}

TriangulatedSurface10::~TriangulatedSurface10()
{
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
  printf("TriangulatedSurface10::volume:\n");
  printf("NOTE: approximate values are used for the interpolated normals\n");

  TriangulatedSurfaceIntegrator tsi(*this);
  RefSCVector norm(_vol->dimension());

  double volume = 0.0;
  RefVertex current = tsi.current();
  for (tsi = 0; tsi; tsi++) {
      tsi.normal(norm);
      volume += tsi.w() * norm[0] * current->point()[0];
    }
  
  // If the volume is negative, then the surface gradients were
  // opposite in sign to the direction assumed.  Flip the sign
  // to fix.
  return fabs(volume);
}

Edge*
TriangulatedSurface10::newEdge(const RefVertex& v0, const RefVertex& v1) const
{
  return new Edge4(v0,v1,_vol,_isovalue);
}

Triangle*
TriangulatedSurface10::newTriangle(const RefEdge& e0,
                                   const RefEdge& e1,
                                   const RefEdge& e2,
                                   int orientation) const
{
  // by construction the edges should be edge4's
  Edge4* e40 = (Edge4*) e0.pointer();
  Edge4* e41 = (Edge4*) e1.pointer();
  Edge4* e42 = (Edge4*) e2.pointer();
  return new Triangle10(e40,e41,e42,_vol,_isovalue,orientation);
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
      RefSCDimension n = _ts->vertex(0)->point().dim();
      RefSCVector p(n);
      RefSCVector g(n);
      _current = new Vertex(p,g);
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

TriangulatedSurfaceIntegrator::
  operator int()
{
  if (_itri < 0 || _itri >= _ts->ntriangle()) return 0;

  TriangleIntegrator* i = _ts->integrator(_itri).pointer();
  _s = i->s(_irs);
  _r = i->r(_irs);
  _weight = i->w(_irs);
  _surface_element = _ts->triangle(_itri)->interpolate(_r,_s,_current);

  //printf("%3d: r=%5.3f s=%5.3f w=%8.5f dA=%8.5f ",
  //       _itri, _r, _s, _weight, _surface_element);
  //_current->print(stdout);

  return (int) 1;
}

void
TriangulatedSurfaceIntegrator::normal(const RefSCVector&n) const
{
  _ts->triangle(_itri)->normal(_r,_s,n);
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

/////////////////////////////////////////////////////////////////////////
// TriangulatedImplicitSurface

#define CLASSNAME TriangulatedImplicitSurface
#define PARENTS public DescribedClass
#define HAVE_KEYVAL_CTOR
//#include <util/state/statei.h>
#include <util/class/classi.h>
void *
TriangulatedImplicitSurface::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

TriangulatedImplicitSurface::
TriangulatedImplicitSurface(const RefKeyVal&keyval)
{
  vol_ = keyval->describedclassvalue("volume");
  if (keyval->error() != KeyVal::OK) {
      fprintf(stderr,"TriangulatedImplicitSurface(const RefKeyVal&keyval): "
              "requires \"volume\"\n");
      abort();
    }

  surf_ = keyval->describedclassvalue("surface");
  if (keyval->error() != KeyVal::OK) surf_ = new TriangulatedSurface;

  isovalue_ = keyval->doublevalue("value");
  if (keyval->error() != KeyVal::OK) isovalue_ = 0.0;

  remove_short_edges_ = keyval->booleanvalue("remove_short_edges");
  if (keyval->error() != KeyVal::OK) remove_short_edges_ = 1;

  remove_slender_triangles_ = keyval->booleanvalue("remove_slender_triangles");
  if (keyval->error() != KeyVal::OK) remove_slender_triangles_ = 0;

  short_edge_factor_ = keyval->doublevalue("short_edge_factor");
  if (keyval->error() != KeyVal::OK) short_edge_factor_ = 0.3;

  slender_triangle_factor_ = keyval->doublevalue("slender_triangle_factor");
  if (keyval->error() != KeyVal::OK) slender_triangle_factor_ = 0.3;

  resolution_ = keyval->doublevalue("resolution");
  if (keyval->error() != KeyVal::OK) resolution_ = 1.0;

  init();
}

void
TriangulatedImplicitSurface::init()
{
  ImplicitSurfacePolygonizer isogen(vol_);
  isogen.set_resolution(resolution_);

  isogen.isosurface(isovalue_,*surf_.pointer());
  surf_->fix_orientation();
  if (remove_short_edges_) {
      surf_->remove_short_edges(short_edge_factor_*resolution_);
      surf_->fix_orientation();
    }
  if (remove_slender_triangles_) {
      surf_->remove_slender_triangles(slender_triangle_factor_);
      surf_->fix_orientation();
    }
}

TriangulatedImplicitSurface::~TriangulatedImplicitSurface()
{
}
