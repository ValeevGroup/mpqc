
#ifndef _util_render_polylines_h
#define _util_render_polylines_h

#include <util/keyval/keyval.h>
#include <util/render/object.h>

class RenderedPolylines: public RenderedObject {
#   define CLASSNAME RenderedPolylines
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  public:
    enum Coloring { None, Vertex };
  private:
    RenderedPolylines::Coloring coloring_;
    int nvertex_;
    int npolyline_;
    int *nvertex_in_polyline_;
    int **polylines_;
    double **vertices_;
    double **vertex_rgb_;
  public:
    RenderedPolylines();
    RenderedPolylines(const RefKeyVal&);
    ~RenderedPolylines();

    void initialize(int nvertex, int npolylines,
                    RenderedPolylines::Coloring c = RenderedPolylines::None);
    int npolyline() { return npolyline_; }
    int nvertex() { return nvertex_; }
    int nvertex_in_polyline(int i) const { return nvertex_in_polyline_[i]; }
    double vertex(int i, int j) const { return vertices_[i][j]; }
    double vertex_rgb(int i, int j) const { return vertex_rgb_[i][j]; }
    int polyline(int i, int j) const { return polylines_[i][j]; }
    int have_vertex_rgb() const { return coloring_ == Vertex; }

    void set_vertex(int, double x, double y, double z);
    void set_vertex_rgb(int, double r, double g, double b);
    void set_polyline(int i, int v1, int v2);
    void set_polyline(int i, int v1, int v2, int v3);
    void set_polyline(int i, int v1, int v2, int v3, int v4);

    void render(const RefRender&);
};
DescribedClass_REF_dec(RenderedPolylines);

#endif
