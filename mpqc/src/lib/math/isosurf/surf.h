//
// surf.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifndef _math_isosurf_surf_h
#define _math_isosurf_surf_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/isosurf/triangle.h>
#include <math/isosurf/volume.h>
#include <util/render/render.h>

class TriangulatedSurface: public DescribedClass {
#   define CLASSNAME TriangulatedSurface
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
//#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    int _verbose;
    int _debug;

    int _completed_surface;

    // sets of objects that make up the surface
    AVLSet<RefVertex> _vertices;
    AVLSet<RefEdge> _edges;
    AVLSet<RefTriangle> _triangles;

    // map objects to an integer index
    AVLMap<RefVertex,int> _vertex_to_index;
    AVLMap<RefEdge,int> _edge_to_index;
    AVLMap<RefTriangle,int> _triangle_to_index;

    // map integer indices to an object
    Array<RefVertex> _index_to_vertex;
    Array<RefEdge> _index_to_edge;
    Array<RefTriangle> _index_to_triangle;

    // mappings between array element numbers
    int** _triangle_vertex;
    int** _triangle_edge;
    int** _edge_vertex;

    // values for each of the vertices
    int _have_values;
    Arraydouble _values;

    // what to use to integrate over the surface, by default
    RefTriangleIntegrator _integrator;
    // other integrators, in terms of time & accuracy:
    // _fast_integrator <= _integrator <= _accurate_interator
    RefTriangleIntegrator _fast_integrator;
    RefTriangleIntegrator _accurate_integrator;

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
    AVLMap<RefVertex,AVLSet<RefEdge> > _tmp_edges;
  public:
    TriangulatedSurface();
    TriangulatedSurface(const RefKeyVal&);
    virtual ~TriangulatedSurface();

    // control printing
    int verbose() const { return _verbose; }
    void verbose(int v) { _verbose = v; }

    // set up an integrator
    void set_integrator(const RefTriangleIntegrator&);
    void set_fast_integrator(const RefTriangleIntegrator&);
    void set_accurate_integrator(const RefTriangleIntegrator&);
    virtual RefTriangleIntegrator integrator(int itri);
    virtual RefTriangleIntegrator fast_integrator(int itri);
    virtual RefTriangleIntegrator accurate_integrator(int itri);

    // construct the surface
    void add_triangle(const RefVertex&,
                      const RefVertex&,
                      const RefVertex&);
    RefEdge find_edge(const RefVertex&, const RefVertex&);
    virtual void complete_surface();

    // clean up the surface
    virtual void remove_short_edges(double cutoff_length = 1.0e-6,
                                    const RefVolume &vol=0, double isoval=0.0);
    virtual void remove_slender_triangles(
                                    int remove_slender, double height_cutoff,
                                    int remove_small, double area_cutoff,
                                    const RefVolume &vol=0, double isoval=0.0);
    virtual void fix_orientation();
    virtual void clear();

    // get information from the object sets
    int nvertex() const { return _vertices.length(); };
    RefVertex vertex(int i) const { return _index_to_vertex[i]; };
    int nedge() const { return _edges.length(); };
    RefEdge edge(int i) const { return _index_to_edge[i]; };
    int ntriangle() const { return _triangles.length(); };
    RefTriangle triangle(int i) const { return _index_to_triangle[i]; }

    // information from the index mappings
    int triangle_vertex(int i,int j) const { return _triangle_vertex[i][j]; };
    int triangle_edge(int i,int j) const { return _triangle_edge[i][j]; };
    int edge_vertex(int i,int j) const { return _edge_vertex[i][j]; };

    // associate values with vertices
    //void compute_colors(Volume&);
    void compute_values(RefVolume&);

    // properties of the surface
    virtual double flat_area(); // use flat triangles
    virtual double flat_volume(); // use flat triangles
    virtual double area();
    virtual double volume();

    // output of the surface
    virtual void print(ostream&o=cout) const;
    virtual void print_vertices_and_triangles(ostream&o=cout) const;
    virtual void print_geomview_format(ostream&o=cout) const;
    virtual void render(const RefRender &render);

    // print information about the topology
    void topology_info(ostream&o=cout);
    void topology_info(int nvertex, int nedge, int ntri, ostream&o=cout);
};
DescribedClass_REF_dec(TriangulatedSurface);

class TriangulatedSurfaceIntegrator {
  private:
    RefTriangulatedSurface _ts;
    int _itri;
    int _irs;
    double _r;
    double _s;
    double _weight;
    double _surface_element;
    RefVertex _current;
    SCVector3 _dA;
    RefTriangleIntegrator (TriangulatedSurface::*_integrator)(int itri);
    RefMessageGrp _grp;
  public:
    TriangulatedSurfaceIntegrator();
    // the surface cannot be changed until this is destroyed
    TriangulatedSurfaceIntegrator(const RefTriangulatedSurface&);
    ~TriangulatedSurfaceIntegrator();
    // Objects initialized by these operators are not automatically
    // updated.  This must be done with the update member.
    // The _grp is not copied.
    void operator = (const TriangulatedSurfaceIntegrator&);
    TriangulatedSurfaceIntegrator(const TriangulatedSurfaceIntegrator&i) {
        operator = (i);
      }
    // Return the number of integration points.
    int n();
    // Assign the surface.  Don't do this while iterating.
    void set_surface(const RefTriangulatedSurface&);
    // returns the number of the vertex in the current triangle
    int vertex_number(int i);
    inline double r() const { return _r; }
    inline double s() const { return _s; }
    inline double w() const { return _weight*_surface_element; }
    double surface_element() const { return _surface_element; }
    double weight() const { return _weight; }
    const SCVector3& dA() const { return _dA; }
    RefVertex current();
    // Tests to see if this point is valid, if it is then
    // _r, _s, etc are computed and 1 is returned.
    int update();
    // This can be used to loop through unique pairs of points.
    // The argument should be a TriangulatedSurfaceIntegrator for
    // the same surface as this.
    int operator < (TriangulatedSurfaceIntegrator&i) {
        update();
        return _itri<i._itri?1:(_itri>i._itri?0:(_irs<i._irs?1:0));
      }
    // Goes to the next point.  Does not update.
    void operator++();
    inline void operator++(int) { operator++(); }
    // setting TSI = i sets TSI to begin at the triangle i
    int operator = (int);
    int itri() const { return _itri; }
    int irs() const { return _irs; }
    // the number of points in the current triangle
    int n_in_tri() const { return (_ts.pointer()->*_integrator)(_itri)->n(); }
    void distribute(const RefMessageGrp &);
    void use_fast_integrator();
    void use_accurate_integrator();
    void use_default_integrator();
};

class TriangulatedImplicitSurface: public TriangulatedSurface {
#   define CLASSNAME TriangulatedImplicitSurface
#   define HAVE_KEYVAL_CTOR
//#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    // The surface is defined as an isosurface of the volume vol_.
    RefVolume vol_;
    double isovalue_;

    int fix_orientation_;
    int remove_short_edges_;
    double short_edge_factor_;
    int remove_slender_triangles_;
    double slender_triangle_factor_;
    int remove_small_triangles_;
    double small_triangle_factor_;
    double resolution_;

    int order_;
  public:
    TriangulatedImplicitSurface(const RefKeyVal&);
    ~TriangulatedImplicitSurface();

    RefVolume volume() const { return vol_; }
    double isovalue() const { return isovalue_; }

    void init();
};
DescribedClass_REF_dec(TriangulatedImplicitSurface);

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
