
#include <stdio.h>
#include <stdlib.h>

#include <util/render/render.h>
#include <util/render/object.h>
#include <util/render/polygons.h>
#include <util/render/color.h>

#define CLASSNAME RenderedPolygons
#define HAVE_KEYVAL_CTOR
#define PARENTS public RenderedObject
#include <util/class/classi.h>
void *
RenderedPolygons::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = RenderedObject::_castdown(cd);
  return do_castdowns(casts,cd);
}

RenderedPolygons::RenderedPolygons()
{
  nvertex_ = 0;
  nface_ = 0;
  vertices_ = 0;
  vertex_rgb_ = 0;
  faces_ = 0;
  nvertex_in_face_ = 0;
}

RenderedPolygons::RenderedPolygons(const RefKeyVal& keyval):
  RenderedObject(keyval)
{
  int nvertex = keyval->count("vertices");
  int nface = keyval->count("faces");
  Coloring coloring;
  if (keyval->count("vertex_color_list")) {
      coloring = Vertex;
    }
  initialize(nvertex, nface, coloring);

  int i;
  for (i=0; i<nvertex; i++) {
      set_vertex(i,
                 keyval->doublevalue("vertices", i, 0),
                 keyval->doublevalue("vertices", i, 1),
                 keyval->doublevalue("vertices", i, 2));
    }

  if (coloring == Vertex) {
      for (i=0; i<nvertex; i++) {
          set_vertex_rgb(i,
                         keyval->doublevalue("vertex_color_list", i, 0),
                         keyval->doublevalue("vertex_color_list", i, 1),
                         keyval->doublevalue("vertex_color_list", i, 2));
        }
    }

  for (i=0; i<nface; i++) {
      nvertex_in_face_[i] = keyval->count("faces", i);
      faces_[i] = new int[nvertex_in_face_[i]];
      for (int j=0; j<nvertex_in_face_[i]; j++) {
          faces_[i][j] = keyval->intvalue("faces", i, j);
        }
    }
}

RenderedPolygons::~RenderedPolygons()
{
  if (vertices_ && vertices_[0]) delete[] vertices_[0];
  if (vertices_) delete[] vertices_;
  if (vertex_rgb_ && vertex_rgb_[0]) delete[] vertex_rgb_[0];
  if (vertex_rgb_) delete[] vertex_rgb_;
  if (faces_) {
      for (int i=0; i<nface_; i++) {
          if (faces_[i]) delete[] faces_[i];
        }
      delete[] faces_;
    }
  if (nvertex_in_face_) {
      delete[] nvertex_in_face_;
    }
}

void
RenderedPolygons::render(const RefRender& render)
{
  render->polygons(this);
}

void
RenderedPolygons::initialize(int nvertex, int nface,
                             RenderedPolygons::Coloring coloring)
{
  coloring_ = coloring;
  nvertex_ = nvertex;
  nface_ = nface;
  
  vertices_ = new double*[nvertex];
  double* tmp = vertices_[0] = new double[3*nvertex];
  int i;
  for (i=1; i<nvertex; i++) {
      tmp += 3;
      vertices_[i] = tmp;
    }

  if (coloring == Vertex) {
      vertex_rgb_ = new double*[nvertex];
      double*tmp = vertex_rgb_[0] = new double[3*nvertex];
      for (int i=1; i<nvertex; i++) {
          tmp += 3;
          vertex_rgb_[i] = tmp;
        }
    }
  else {
      vertex_rgb_ = 0;
    }

  faces_ = new int*[nface];
  nvertex_in_face_ = new int[nface];
  for (i=0; i<nface; i++) {
      faces_[i] = 0;
      nvertex_in_face_[i] = 0;
    }
}

void
RenderedPolygons::set_vertex(int i, double x, double y, double z)
{
  vertices_[i][0] = x;
  vertices_[i][1] = y;
  vertices_[i][2] = z;
}

void
RenderedPolygons::set_vertex_rgb(int i, double r, double g, double b)
{
  vertex_rgb_[i][0] = r;
  vertex_rgb_[i][1] = g;
  vertex_rgb_[i][2] = b;
}

void
RenderedPolygons::set_face(int iface, int v1, int v2, int v3)
{
  if (faces_[iface]) delete[] faces_[iface];
  faces_[iface] = new int[3];
  faces_[iface][0] = v1;
  faces_[iface][1] = v2;
  faces_[iface][2] = v3;

  nvertex_in_face_[iface] = 3;
}

void
RenderedPolygons::set_face(int iface, int v1, int v2, int v3, int v4)
{
  if (faces_[iface]) delete[] faces_[iface];
  faces_[iface] = new int[4];
  faces_[iface][0] = v1;
  faces_[iface][1] = v2;
  faces_[iface][2] = v3;
  faces_[iface][3] = v4;

  nvertex_in_face_[iface] = 4;
}
