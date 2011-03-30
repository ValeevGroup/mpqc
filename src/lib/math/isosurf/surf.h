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

#ifdef HAVE_CONFIG_H
#include <scconfig.h>
#endif

#include <map>
#include <set>
#include <vector>

#include <math/isosurf/triangle.h>
#include <math/isosurf/volume.h>
#include <util/render/render.h>

namespace sc {

template <class C, class I>
inline void
erase_elements_by_value(C &container, I begin, I end)
{
  for (I i=begin; i!=end; i++) {
      container.erase(*i);
    }
}

class TriangulatedSurface: public DescribedClass {
  protected:
    int _verbose;
    int _debug;

    int _completed_surface;

    // sets of objects that make up the surface
    std::set<Ref<Vertex> > _vertices;
    std::set<Ref<Edge> > _edges;
    std::set<Ref<Triangle> > _triangles;

    // map objects to an integer index
    std::map<Ref<Vertex>,int> _vertex_to_index;
    std::map<Ref<Edge>,int> _edge_to_index;
    std::map<Ref<Triangle>,int> _triangle_to_index;

    // map integer indices to an object
    std::vector<Ref<Vertex> > _index_to_vertex;
    std::vector<Ref<Edge> > _index_to_edge;
    std::vector<Ref<Triangle> > _index_to_triangle;

    // mappings between array element numbers
    int** _triangle_vertex;
    int** _triangle_edge;
    int** _edge_vertex;

    // values for each of the vertices
    int _have_values;
    std::vector<double> _values;

    // what to use to integrate over the surface, by default
    Ref<TriangleIntegrator> _integrator;
    // other integrators, in terms of time & accuracy:
    // _fast_integrator <= _integrator <= _accurate_interator
    Ref<TriangleIntegrator> _fast_integrator;
    Ref<TriangleIntegrator> _accurate_integrator;

    void clear_int_arrays();

    void complete_ref_arrays();
    void complete_int_arrays();

    void recompute_index_maps();

    void add_triangle(const Ref<Triangle>&);
    void add_vertex(const Ref<Vertex>&);
    void add_edge(const Ref<Edge>&);

    // these members must be used to allocate new triangles and edges
    // since specializations of TriangulatedSurface might need to
    // override these to produce triangles and edges with interpolation
    // data.
    virtual Triangle* newTriangle(const Ref<Edge>&,
                                  const Ref<Edge>&,
                                  const Ref<Edge>&,
                                  int orientation) const;
    virtual Edge* newEdge(const Ref<Vertex>&,const Ref<Vertex>&) const;

    // this map of edges to vertices is used to construct the surface
    std::map<Ref<Vertex>,std::set<Ref<Edge> > > _tmp_edges;
  public:
    TriangulatedSurface();
    TriangulatedSurface(const Ref<KeyVal>&);
    virtual ~TriangulatedSurface();

    // control printing
    int verbose() const { return _verbose; }
    void verbose(int v) { _verbose = v; }

    // set up an integrator
    void set_integrator(const Ref<TriangleIntegrator>&);
    void set_fast_integrator(const Ref<TriangleIntegrator>&);
    void set_accurate_integrator(const Ref<TriangleIntegrator>&);
    virtual Ref<TriangleIntegrator> integrator(int itri);
    virtual Ref<TriangleIntegrator> fast_integrator(int itri);
    virtual Ref<TriangleIntegrator> accurate_integrator(int itri);

    // construct the surface
    void add_triangle(const Ref<Vertex>&,
                      const Ref<Vertex>&,
                      const Ref<Vertex>&);
    Ref<Edge> find_edge(const Ref<Vertex>&, const Ref<Vertex>&);
    virtual void complete_surface();

    // clean up the surface
    virtual void remove_short_edges(double cutoff_length = 1.0e-6,
                                    const Ref<Volume> &vol=0, double isoval=0.0);
    virtual void remove_slender_triangles(
                                    int remove_slender, double height_cutoff,
                                    int remove_small, double area_cutoff,
                                    const Ref<Volume> &vol=0, double isoval=0.0);
    virtual void fix_orientation();
    virtual void clear();

    // get information from the object sets
    int nvertex() const { return _vertices.size(); };
    Ref<Vertex> vertex(int i) const { return _index_to_vertex[i]; };
    int vertex_index(const Ref<Vertex> &o) {
      std::map<Ref<Vertex>,int>::iterator i = _vertex_to_index.find(o);
      if (i != _vertex_to_index.end()) return i->second;
      return -1;
    }
    int nedge() const { return _edges.size(); };
    Ref<Edge> edge(int i) const { return _index_to_edge[i]; };
    int edge_index(const Ref<Edge> &o) {
      std::map<Ref<Edge>,int>::iterator i = _edge_to_index.find(o);
      if (i != _edge_to_index.end()) return i->second;
      return -1;
    }
    int ntriangle() const { return _triangles.size(); };
    Ref<Triangle> triangle(int i) const { return _index_to_triangle[i]; }
    int triangle_index(const Ref<Triangle> &o) {
      std::map<Ref<Triangle>,int>::iterator i = _triangle_to_index.find(o);
      if (i != _triangle_to_index.end()) return i->second;
      return -1;
    }

    // information from the index mappings
    int triangle_vertex(int i,int j) const { return _triangle_vertex[i][j]; };
    int triangle_edge(int i,int j) const { return _triangle_edge[i][j]; };
    int edge_vertex(int i,int j) const { return _edge_vertex[i][j]; };

    // associate values with vertices
    //void compute_colors(Volume&);
    void compute_values(Ref<Volume>&);

    // properties of the surface
    virtual double flat_area(); // use flat triangles
    virtual double flat_volume(); // use flat triangles
    virtual double area();
    virtual double volume();

    // output of the surface
    virtual void print(std::ostream&o=ExEnv::out0()) const;
    virtual void print_vertices_and_triangles(std::ostream&o=ExEnv::out0()) const;
    virtual void print_geomview_format(std::ostream&o=ExEnv::out0()) const;
    virtual void render(const Ref<Render> &render);

    // print information about the topology
    void topology_info(std::ostream&o=ExEnv::out0());
    void topology_info(int nvertex, int nedge, int ntri, std::ostream&o=ExEnv::out0());
};


class TriangulatedSurfaceIntegrator {
  private:
    Ref<TriangulatedSurface> _ts;
    int _itri;
    int _irs;
    double _r;
    double _s;
    double _weight;
    double _surface_element;
    Ref<Vertex> _current;
    SCVector3 _dA;
    Ref<TriangleIntegrator> (TriangulatedSurface::*_integrator)(int itri);
    Ref<MessageGrp> _grp;
  public:
    TriangulatedSurfaceIntegrator();
    // the surface cannot be changed until this is destroyed
    TriangulatedSurfaceIntegrator(const Ref<TriangulatedSurface>&);
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
    void set_surface(const Ref<TriangulatedSurface>&);
    // returns the number of the vertex in the current triangle
    int vertex_number(int i);
    inline double r() const { return _r; }
    inline double s() const { return _s; }
    inline double w() const { return _weight*_surface_element; }
    double surface_element() const { return _surface_element; }
    double weight() const { return _weight; }
    const SCVector3& dA() const { return _dA; }
    Ref<Vertex> current();
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
    void distribute(const Ref<MessageGrp> &);
    void use_fast_integrator();
    void use_accurate_integrator();
    void use_default_integrator();
};

class TriangulatedImplicitSurface: public TriangulatedSurface {
  private:
    // The surface is defined as an isosurface of the volume vol_.
    Ref<Volume> vol_;
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

    int inited_;
  public:
    TriangulatedImplicitSurface(const Ref<KeyVal>&);
    ~TriangulatedImplicitSurface();

    Ref<Volume> volume_object() const { return vol_; }
    double isovalue() const { return isovalue_; }

    void init();
    int inited() const { return inited_; }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
