
#ifndef _math_isosurf_surf_h
#define _math_isosurf_surf_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/isosurf/triangle.h>

class TriangulatedSurface {
  protected:
    int _completed_surface;

    // sets of objects that make up the surface
    ArraysetRefVertex _vertices;
    ArraysetRefEdge _edges;
    ArraysetRefTriangle _triangles;

    // mappings between array element numbers
    int** _triangle_vertex;
    int** _triangle_edge;
    int** _edge_vertex;

    // values for each of the vertices
    int _have_values;
    Arraydouble _values;

    // what to use to integrate over the surface
    TriangleIntegrator* _integrator;
  public:
    TriangulatedSurface();
    virtual ~TriangulatedSurface();

    // set up an integrator
    void set_integrator(TriangleIntegrator*); // arg is deleted by DTOR or set
    virtual TriangleIntegrator* integrator(int itri);

    // construct the surface
    void add_triangle(const RefTriangle&);
    RefEdge find_edge(const RefVertex&, const RefVertex&);
    virtual void complete_surface();

    // clean up the surface
    virtual void remove_short_edges(double cutoff_length = 1.0e-6);
    virtual void clear();

    // get information from the object sets
    inline int nvertex() { return _vertices.length(); };
    inline RefVertex vertex(int i) { return _vertices[i]; };
    inline int nedge() { return _edges.length(); };
    inline RefEdge edge(int i) { return _edges[i]; };
    inline int ntriangle() { return _triangles.length(); };
    inline RefTriangle triangle(int i) { return _triangles[i]; };

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
};

class TriangulatedSurface10: public TriangulatedSurface {
  private:
    RefVolume _vol;
    double _isovalue;
    ArrayRefEdge4 _edges4;
    ArrayRefTriangle10 _triangles10;
  public:
    TriangulatedSurface10(RefVolume&vol,double isovalue);
    ~TriangulatedSurface10();
    void remove_short_edges(double cutoff_length = 1.0e-6);
    void complete_surface();
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
    inline double r() { return _r; }
    inline double s() { return _s; }
    inline double w() { return _weight*_surface_element; }
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

#endif
