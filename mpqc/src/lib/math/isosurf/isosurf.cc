
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
  _volume(vol)
{
}

ImplicitSurfacePolygonizer::~ImplicitSurfacePolygonizer()
{
}

ImplicitSurfacePolygonizer* ImplicitSurfacePolygonizer::current = 0;

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

//   fprintf(stderr,"bounding box is (%f, %f, %f) (%f, %f, %f)\n",
//           p0.get_element(0),p0.get_element(1),p0.get_element(2),
//           p1.get_element(0),p1.get_element(1),p1.get_element(2));
//   fprintf(stderr,"midpoint is (%f, %f, %f)\n",
//           midpoint.get_element(0),
//           midpoint.get_element(1),
//           midpoint.get_element(2));

  // Initialize the static pointer to this, so the C polygonizer can find us.
  current = this;
  _surf = &surf;
  _value = value;
  // Find the polygons.
  char *msg = polygonize(value_of_current, _resolution, bounds,
                         midpoint(0), midpoint(1), midpoint(2),
                         add_triangle_to_current, NOTET);
  current = 0;
  _surf = 0;
  if (msg) {
      fprintf(stderr, "ImplicitSurfacePolygonizer::isosurface: failed: %s\n",
              msg);
      abort();
    }

  // Clean up temporaries.
  _tmp_vertices.clear();

  // finish the surface
  surf.complete_surface();
}

double
ImplicitSurfacePolygonizer::value_of_current(double x,double y,double z)
{
  RefSCVector v(current->_volume->dimension());
  v(0) = x; v(1) = y; v(2) = z;
  current->_volume->set_x(v);
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

  RefEdge e0 = current->_surf->find_edge(v1,v2);
  RefEdge e1 = current->_surf->find_edge(v2,v3);
  RefEdge e2 = current->_surf->find_edge(v3,v1);

  if (e0.null()) {
      e0 = new Edge(v1,v2);
    }
  if (e1.null()) {
      e1 = new Edge(v2,v3);
    }
  if (e2.null()) {
      e2 = new Edge(v3,v1);
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
