
#ifndef _math_isosurf_surf_h
#define _math_isosurf_surf_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/isosurf/triangle.h>
#include <math/isosurf/vertexAVLSet.h>
#include <math/isosurf/edgeAVLSet.h>
#include <math/isosurf/triAVLSet.h>
#include <util/container/pixintRAVLMap.h>
#include <util/container/intpixRAVLMap.h>
#include <math/isosurf/edgeRAVLMap.h>
#include <math/isosurf/volume.h>

class TriangulatedSurface: public DescribedClass {
#   define CLASSNAME TriangulatedSurface
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
//#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    int _completed_surface;

    // sets of objects that make up the surface
    RefVertexAVLSet _vertices;
    RefEdgeAVLSet _edges;
    RefTriangleAVLSet _triangles;

    // map pixes to an integer index
    PixintRAVLMap _vertex_to_index;
    PixintRAVLMap _edge_to_index;
    PixintRAVLMap _triangle_to_index;

    // map integer indices to a pix
    intPixRAVLMap _index_to_vertex;
    intPixRAVLMap _index_to_edge;
    intPixRAVLMap _index_to_triangle;

    // mappings between array element numbers
    int** _triangle_vertex;
    int** _triangle_edge;
    int** _edge_vertex;

    // values for each of the vertices
    int _have_values;
    Arraydouble _values;

    // what to use to integrate over the surface
    RefTriangleIntegrator _integrator;

    void clear_int_arrays();

    void complete_ref_arrays();
    void complete_int_arrays();

    void recompute_index_maps();

    void add_triangle(const RefTriangle&);
    void add_vertex(const RefVertex&);
    void add_edge(const RefEdge&);

    // these members must be used to allocate new triangles and edges
    // since specializations of TriangulatedSurface might need to
    // override these to produce triangles and edges with interpolation
    // data.
    virtual Triangle* newTriangle(const RefEdge&,
                                  const RefEdge&,
                                  const RefEdge&,
                                  int orientation) const;
    virtual Edge* newEdge(const RefVertex&,const RefVertex&) const;

    // this map of edges to vertices is used to construct the surface
    RefVertexRefEdgeAVLSetRAVLMap _tmp_edges;
  public:
    TriangulatedSurface();
    TriangulatedSurface(const RefKeyVal&);
    virtual ~TriangulatedSurface();

    // set up an integrator
    void set_integrator(const RefTriangleIntegrator&);
    virtual RefTriangleIntegrator integrator(int itri);

    // construct the surface
    void add_triangle(const RefVertex&,
                      const RefVertex&,
                      const RefVertex&);
    RefEdge find_edge(const RefVertex&, const RefVertex&);
    virtual void complete_surface();

    // clean up the surface
    virtual void remove_short_edges(double cutoff_length = 1.0e-6);
    virtual void remove_slender_triangles(double heigth_cutoff = 1.0e-6);
    virtual void fix_orientation();
    virtual void clear();

    // get information from the object sets
    inline int nvertex() { return _vertices.length(); };
    inline RefVertex vertex(int i) { return _vertices(_index_to_vertex[i]); };
    inline int nedge() { return _edges.length(); };
    inline RefEdge edge(int i) { return _edges(_index_to_edge[i]); };
    inline int ntriangle() { return _triangles.length(); };
    inline RefTriangle triangle(int i) {
        return _triangles(_index_to_triangle[i]);
      }

    // information from the index mappings
    inline int triangle_vertex(int i,int j) { return _triangle_vertex[i][j]; };
    inline int triangle_edge(int i,int j) { return _triangle_edge[i][j]; };
    inline int edge_vertex(int i,int j) { return _edge_vertex[i][j]; };

    // associate values with vertices
    //void compute_colors(Volume&);
    void compute_values(RefVolume&);

    // properties of the surface
    virtual double area();
    virtual double volume();

    // output of the surface
    virtual void print(FILE* = stdout);
    virtual void print_vertices_and_triangles(FILE* = stdout);
    virtual void print_geomview_format(FILE*fp = stdout);

    // print information about the topology
    void topology_info(FILE*f=stdout);
    void topology_info(int nvertex, int nedge, int ntri, FILE*f=stdout);
};
DescribedClass_REF_dec(TriangulatedSurface);

class TriangulatedSurface10: public TriangulatedSurface {
#   define CLASSNAME TriangulatedSurface10
#   define HAVE_KEYVAL_CTOR
//#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    RefVolume _vol;
    double _isovalue;
    ArrayRefEdge4 _edges4;
    ArrayRefTriangle10 _triangles10;
  protected:
    Triangle* newTriangle(const RefEdge&,
                          const RefEdge&,
                          const RefEdge&,
                          int orientation) const;
    Edge* newEdge(const RefVertex&,const RefVertex&) const;
  public:
    TriangulatedSurface10(const RefVolume&vol,double isovalue);
    TriangulatedSurface10(const RefKeyVal& keyval);
    ~TriangulatedSurface10();
    double area();
    double volume();
};

class TriangulatedSurfaceIntegrator {
  private:
    TriangulatedSurface* _ts;
    int _itri;
    int _irs;
    double _r;
    double _s;
    double _weight;
    double _surface_element;
    RefVertex _current;
  public:
    // the surface cannot be changed until this is destroyed
    TriangulatedSurfaceIntegrator(TriangulatedSurface&);
    ~TriangulatedSurfaceIntegrator();
    // returns the number of the vertex in the current triangle
    int vertex_number(int i);
    inline double r() const { return _r; }
    inline double s() const { return _s; }
    inline double w() const { return _weight*_surface_element; }
    void normal(const RefSCVector&) const;
    RefVertex current();
    // Tests to see if this point is value, if it is then
    // _r, _s, etc are computed.
    // NOTE: this is the only member that can cause _r, _s, etc
    // to be computed.
    operator int();
    // go to the next point
    void operator++();
    inline void operator++(int) { operator++(); }
    // setting TSI = i sets TSI to begin at the triangle i
    int operator = (int);
};

class TriangulatedImplicitSurface: public DescribedClass {
#   define CLASSNAME TriangulatedImplicitSurface
#   define HAVE_KEYVAL_CTOR
//#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    // The surface is defined as an isosurface of the volume vol_.
    RefVolume vol_;
    double isovalue_;
    RefTriangulatedSurface surf_;

    int remove_short_edges_;
    double short_edge_factor_;
    int remove_slender_triangles_;
    double slender_triangle_factor_;
    double resolution_;

    void init();
  public:
    TriangulatedImplicitSurface(const RefKeyVal&);
    ~TriangulatedImplicitSurface();

    int nvertex() { return surf_->nvertex(); }
    RefVertex vertex(int i) { return surf_->vertex(i); }
    int nedge() { return surf_->nedge(); }
    RefEdge edge(int i) { return surf_->edge(i); }
    int ntriangle() { return surf_->ntriangle(); }
    RefTriangle triangle(int i) { return surf_->triangle(i); }

    int triangle_vertex(int i,int j) { return surf_->triangle_vertex(i,j); }
    int triangle_edge(int i,int j) { return surf_->triangle_edge(i,j); }
    int edge_vertex(int i,int j) { return surf_->edge_vertex(i,j); }
};
DescribedClass_REF_dec(TriangulatedImplicitSurface);

#endif
