
#ifndef _util_render_polygons_h
#define _util_render_polygons_h

#include <util/render/object.h>

class RenderedPolygons: public RenderedObject {
#   define CLASSNAME RenderedPolygons
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  private:
    enum Color { None, Vertex /*, Face*/ };
    RenderedPolygons::Color color_;
    int nvertex_;
    int nface_;
    double** vertices_;
    double** vertex_rgb_;
    int** faces_;
    int* nvertex_in_face_;
  protected:
    void render(const RefRender&);
  public:
    RenderedPolygons();
    RenderedPolygons(const RefKeyVal&);
    ~RenderedPolygons();

    void initialize(int nvertex, int nface,
                    RenderedPolygons::Color c = RenderedPolygons::None);
    void set_vertex(int, double x, double y, double z);
    void set_vertex_rgb(int, double r, double g, double b);
    void set_face(int iface, int v1, int v2, int v3);
    void set_face(int iface, int v1, int v2, int v3, int v4);

    int nvertex() const { return nvertex_; }
    int nface() const { return nface_; }
    int nvertex_in_face(int iface) const { return nvertex_in_face_[iface]; }
    double vertex(int i, int j) const { return vertices_[i][j]; }
    int face(int i,int j) const { return faces_[i][j]; }
    double vertex_rgb(int i, int j) const { return vertex_rgb_[i][j]; }
    int have_vertex_rgb() const { return color_ == Vertex; }
};
DescribedClass_REF_dec(RenderedPolygons);

#endif
