
#ifdef __GNUC__
#pragma implementation
#endif

#include <math/isosurf/isosurf.h>
#include <math/isosurf/implicit.h>

////////////////////////////////////////////////////////////////////////////
// IsosurfaceGen members

IsosurfaceGen::IsosurfaceGen():
  _resolution(0.25)
{
}

IsosurfaceGen::~IsosurfaceGen()
{
}

void
IsosurfaceGen::set_resolution(double r)
{
  _resolution = r;
}

////////////////////////////////////////////////////////////////////////////
// ImplicitSurfacePolygonizer members

ImplicitSurfacePolygonizer::ImplicitSurfacePolygonizer(const RefVolume&vol):
  _volume(vol),
  _tmp_edges(RefEdgeAVLSet())
{
}

ImplicitSurfacePolygonizer::~ImplicitSurfacePolygonizer()
{
}

ImplicitSurfacePolygonizer* ImplicitSurfacePolygonizer::current = 0;

static RefSCVector current_x;
void
ImplicitSurfacePolygonizer::isosurface(double value,
                                       TriangulatedSurface& surf)
{
  surf.clear();
  RefSCDimension dim = _volume->dimension();

  // Only 3D volumes are handled.
  if (dim.n() != 3) {
      fprintf(stderr,"ImplicitSurfacePolygonizer::isosurface: "
              " volume dimension is not 3\n");
      abort();
    }

  // Find the bounding box.
  RefSCVector p0(dim);
  RefSCVector p1(dim);
  _volume->boundingbox(value - 0.001, value + 0.001, p0, p1);
  RefSCVector diag = p1 - p0;
  RefSCVector midpoint = 0.5*diag + p0;
  double biggest_width = diag.maxabs();
  int bounds = (int)(0.5*biggest_width/_resolution) + 2;

  // polygonize will find a starting point and do bounds checking
  // from that point.  To make sure the bounding box is big enough
  // its size must be doubled.  Since polygonization is implicit
  // there is no performance penalty.
  bounds *= 2;

  fprintf(stderr,"bounding box is (%f, %f, %f) (%f, %f, %f)\n",
          p0.get_element(0),p0.get_element(1),p0.get_element(2),
          p1.get_element(0),p1.get_element(1),p1.get_element(2));
  fprintf(stderr,"midpoint is (%f, %f, %f)\n",
          midpoint.get_element(0),
          midpoint.get_element(1),
          midpoint.get_element(2));
  fprintf(stderr,"biggest_width = %f, resolution = %f, bounds = %d\n",
          biggest_width, _resolution, bounds);

  // Initialize the static pointer to this, so the C polygonizer can find us.
  current_x = RefSCVector(_volume->dimension());
  current = this;
  _surf = &surf;
  _value = value;
  // Find the polygons.
  char *msg = polygonize(value_of_current, _resolution, bounds,
                         midpoint(0), midpoint(1), midpoint(2),
                         add_triangle_to_current, NOTET);
  current = 0;
  _surf = 0;
  current_x = 0;
  if (msg) {
      fprintf(stderr, "ImplicitSurfacePolygonizer::isosurface: failed: %s\n",
              msg);
      abort();
    }

  // Clean up temporaries.
  _tmp_vertices.clear();
  _tmp_edges.clear();

  // finish the surface
  surf.complete_surface();
}

double
ImplicitSurfacePolygonizer::value_of_current(double x,double y,double z)
{
  current_x(0) = x; current_x(1) = y; current_x(2) = z;
  current->_volume->set_x(current_x);
  return current->_volume->value() - current->_value;
}

int
ImplicitSurfacePolygonizer::add_triangle_to_current(int i1, int i2, int i3,
                                                   VERTICES v)
{
  int i;
  for (i=current->_tmp_vertices.length(); i<v.count; i++) {
      RefSCVector newpoint(current->_volume->dimension());
      newpoint[0] = v.ptr[i].position.x;
      newpoint[1] = v.ptr[i].position.y;
      newpoint[2] = v.ptr[i].position.z;
      current->_volume->set_x(newpoint);
      RefSCVector gradient = current->_volume->gradient();
      current->_tmp_vertices.add(new Vertex(newpoint, gradient));
    }

  RefVertex v1 = current->_tmp_vertices[i1];
  RefVertex v2 = current->_tmp_vertices[i2];
  RefVertex v3 = current->_tmp_vertices[i3];

  // Find this triangles edges if they have already be created
  // for some other triangle.
#if 0
  // this is rather slow
  RefEdge e0 = current->_surf->find_edge(v1,v2);
  RefEdge e1 = current->_surf->find_edge(v2,v3);
  RefEdge e2 = current->_surf->find_edge(v3,v1);
#else
  RefEdge e0, e1, e2;
  RefEdgeAVLSet& v1edges = current->_tmp_edges[v1];
  RefEdgeAVLSet& v2edges = current->_tmp_edges[v2];

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
#endif

  if (e0.null()) {
      e0 = new Edge(v1,v2);
      current->_tmp_edges[v1].add(e0);
      current->_tmp_edges[v2].add(e0);
    }
  if (e1.null()) {
      e1 = new Edge(v2,v3);
      current->_tmp_edges[v2].add(e1);
      current->_tmp_edges[v3].add(e1);
    }
  if (e2.null()) {
      e2 = new Edge(v3,v1);
      current->_tmp_edges[v3].add(e2);
      current->_tmp_edges[v1].add(e2);
    }
  
  int orientation;
  if (e0->vertex(0) == v1) {
      orientation = 0;
    }
  else {
      orientation = 1;
    }
  
  current->_surf->add_triangle(new Triangle(e0,e1,e2,orientation));

  return 1;
}
