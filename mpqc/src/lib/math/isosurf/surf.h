
#ifndef _math_isosurf_surf_h
#define _math_isosurf_surf_h

#ifdef __GNUC__
#pragma interface
#endif

extern "C" {
#include <stdio.h>
}
#include <util/container/ref.h>
#include <util/container/set.h>
#include <math/scmat/matrix.h>
#include <math/isosurf/volume.h>

class Vertex: public RefCount {
  private:
    RefSCVector _point;
    RefSCVector _gradient;
  public:
    Vertex();
    Vertex(RefSCVector&point,RefSCVector&gradient);
    ~Vertex();
    RefSCVector gradient();
    RefSCVector point();
    void set_point(RefSCVector&p);
    operator RefSCVector();
    RefSCDimension dimension();

    void print(FILE*fp = stdout);
};

REF_dec(Vertex);
ARRAY_dec(RefVertex);
SET_dec(RefVertex);
ARRAYSET_dec(RefVertex);
ARRAY_dec(ArraysetRefVertex);

class Edge: public VRefCount {
  private:
    RefVertex _vertices[2];
  public:
    inline Edge(RefVertex p1,RefVertex p2) { _vertices[0]=p1; _vertices[1]=p2; };
    inline RefVertex vertex(int i) { return _vertices[i]; };
    virtual ~Edge();
    virtual double length();
    void add_vertices(SetRefVertex&);
};

REF_dec(Edge);
ARRAY_dec(RefEdge);
SET_dec(RefEdge);
ARRAYSET_dec(RefEdge);
ARRAY_dec(ArraysetRefEdge);

// a 4 vertex edge (used in 10 vertex triangles)
class Edge4: public Edge {
  private:
    RefVertex _interiorvertices[2];
  public:
    Edge4(RefVertex p1,RefVertex p2,RefVolume vol,double isovalue);
    inline RefVertex interior_vertex(int i) { return _interiorvertices[i]; };
    ~Edge4();
    double length();
};

REF_dec(Edge4);
ARRAY_dec(RefEdge4);
SET_dec(RefEdge4);
ARRAYSET_dec(RefEdge4);
ARRAY_dec(ArraysetRefEdge4);

class Triangle: public VRefCount {
  private:
    int _edge_vertex;
    RefEdge _edges[3];
  public:
    Triangle(RefEdge v1, RefEdge v2, RefEdge v3);
    inline RefEdge edge(int i) { return _edges[i]; };
    virtual ~Triangle();
    void add_edges(SetRefEdge&);
    void add_vertices(SetRefVertex&);

    // returns the surface area element
    // 0<=r<=1, 0<=s<=1, 0<=r+s<=1
    // RefVertex is the intepolated vertex (both point and gradient)
    virtual double interpolate(double r,double s,RefVertex&v);

    // returns a vertex in the triangle
    // i = 0 is the (0,0) vertex
    // i = 1 is the (r=0,s=1) vertex
    // i = 2 is the (r=1,s=0) vertex
    RefVertex vertex(int i);

    virtual double area();
};

REF_dec(Triangle);
ARRAY_dec(RefTriangle);
SET_dec(RefTriangle);
ARRAYSET_dec(RefTriangle);

// 10 vertex triangle
class Triangle10: public Triangle {
  private:
    RefVertex _vertices[10];
  public:
    Triangle10(RefEdge4 v1, RefEdge4 v2, RefEdge4 v3,
               RefVolume vol,
               double isovalue);
    virtual ~Triangle10();
    double interpolate(double r,double s,RefVertex&v);
    double area();
};

REF_dec(Triangle10);
ARRAY_dec(RefTriangle10);
SET_dec(RefTriangle10);
ARRAYSET_dec(RefTriangle10);

class TriangleIntegrator {
  private:
    int _n;
    double* _r;
    double* _s;
    double* _w;
  protected:
    void set_r(int i,double r);
    void set_s(int i,double s);
    void set_w(int i,double w);
  public:
    TriangleIntegrator(int n);
    virtual ~TriangleIntegrator();
    inline double w(int i) { return _w[i]; }
    inline double r(int i) { return _r[i]; }
    inline double s(int i) { return _s[i]; }
    inline int n() { return _n; }
};

class GaussTriangleIntegrator: public TriangleIntegrator {
  public:
    GaussTriangleIntegrator(int order);
    ~GaussTriangleIntegrator();
};

class TriangulatedSurface {
  protected:
    int _completed_surface;
    
    ArraysetRefVertex _vertices;
    ArraysetRefEdge _edges;
    ArraysetRefTriangle _triangles;

    int** _triangle_vertex;
    int** _triangle_edge;
    int** _edge_vertex;

    int _have_values;
    Arraydouble _values;

    TriangleIntegrator* _integrator;
  public:
    TriangulatedSurface();
    virtual ~TriangulatedSurface();
    void set_integrator(TriangleIntegrator*); // arg is deleted by TS
    virtual TriangleIntegrator* integrator(int itri);
    inline int nvertex() { return _vertices.length(); };
    inline RefVertex vertex(int i) { return _vertices[i]; };
    inline int nedge() { return _edges.length(); };
    inline RefEdge edge(int i) { return _edges[i]; };
    inline int ntriangle() { return _triangles.length(); };
    inline RefTriangle triangle(int i) { return _triangles[i]; };
    void add_triangle(RefTriangle&);
    void initialize_vertices_triangles(int nvertex, int ntriangle);
    void add_vertex(RefVertex&);
    void add_triangle(RefTriangle&,int vertex0,int vertex1,int vertex2);
    virtual void remove_short_edges(double cutoff_length = 1.0e-6);
    virtual void clear();
    virtual void complete_surface();
    virtual void print(FILE* = stdout);
    virtual void print_vertices_and_triangles(FILE* = stdout);
    //void compute_colors(Volume&);
    void compute_values(RefVolume&);
    virtual double area();
    virtual double volume();
    inline int triangle_vertex(int i,int j) { return _triangle_vertex[i][j]; };
    inline int triangle_edge(int i,int j) { return _triangle_edge[i][j]; };
    inline int edge_vertex(int i,int j) { return _edge_vertex[i][j]; };
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

class UniformLattice {
  protected:
    double* _start;
    double* _incr;
    RefSCDimension _scdim;
    int _ndim;
    int* _dim;
  public:
    UniformLattice(RefSCDimension&);
    // 3D CTOR
    UniformLattice(int dim0,double start0,double incr0,
                   int dim1,double start1,double incr1,
                   int dim2,double start2,double incr2);
    virtual ~UniformLattice();
    void set_lattice_parameters(int dim,int*n,double*start,double*incr);
    inline const int* dim_vector() { return _dim; };
    inline RefSCDimension scdim() { return _scdim; };
    inline int dim(int i) { return _dim[i]; };
    inline int ndim() { return _ndim; };
    inline double start(int i) { return _start[i]; };
    inline double incr(int i) { return _incr[i]; };

    // 3D access functions
    double value(int,int,int);
    RefSCVector interpolate(int,int,int,int,int,int,double);

    // general access functions
    virtual double value(int*) = 0;
    virtual RefSCVector interpolate(int*,int*,double) = 0;

    virtual void gradient(RefSCVector&,RefSCVector&) = 0;
};

class ImplicitUniformLattice: public UniformLattice {
  protected:
    RefVolume _vol;
  public:
    ImplicitUniformLattice(RefVolume& volume,
                           double resolution,
                           double minval,
                           double maxval);
    // 3D CTOR
    ImplicitUniformLattice(RefVolume& volume,
                           int dim0,double start0,double incr0,
                           int dim1,double start1,double incr1,
                           int dim2,double start2,double incr2);
    virtual ~ImplicitUniformLattice();

    // 3D access functions
    virtual double value(int*);
    virtual RefSCVector interpolate(int*,int*,double);

    virtual void gradient(RefSCVector&,RefSCVector&);
};

class StoredUniformLattice: public UniformLattice {
  protected:
    double* _data;
  public:
    StoredUniformLattice(int dim0,double start0,double incr0,
                         int dim1,double start1,double incr1,
                         int dim2,double start2,double incr2);
    ~StoredUniformLattice();
    void set_value(int,int,int,double);
    inline double* data() { return _data; };

    double value(int*);
};

class IsosurfaceGen {
  public:
    IsosurfaceGen();
    virtual ~IsosurfaceGen();
    virtual void isosurface(double value,
                            TriangulatedSurface& surf) = 0;
};

class MCubesIsosurfaceGen: public IsosurfaceGen {
  protected:
    UniformLattice& _lattice;
  public:
    MCubesIsosurfaceGen(UniformLattice&);
    virtual ~MCubesIsosurfaceGen();
    virtual void isosurface(double value,
                            TriangulatedSurface& surf);
};  

#endif
