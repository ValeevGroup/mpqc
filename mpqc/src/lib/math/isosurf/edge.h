
#ifndef _math_isosurf_edge_h
#define _math_isosurf_edge_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/isosurf/vertex.h>

class Edge: public VRefCount {
  private:
    int _order;
    RefVertex *_vertices; // nvertices = _order + 1
  public:
    Edge(const RefVertex &p1,
         const RefVertex &p2)
    {
      _order = 1;
      _vertices = new RefVertex[2];
      _vertices[0]=p1; _vertices[1]=p2;
    }
    Edge(const RefVertex &p1,
         const RefVertex &p2,
         const RefVertex &p3);
    Edge(const RefVertex &p1,
         const RefVertex &p2,
         const RefVertex &p3,
         const RefVertex &p4);
    ~Edge();
    int order() const { return _order; }
    double straight_length();
    // return the endpoints
    RefVertex vertex(int i) const
    {
      return i?_vertices[_order]:_vertices[0];
    }
    // returns endpoints or interior vertex 0 <= i <= order
    RefVertex interior_vertex(int i) const
    {
      return _vertices[i];
    }
    // add the endpoints to the set
    void add_vertices(SetRefVertex&);
    void set_order(int order, const RefVolume&vol,double isovalue);
};

REF_dec(Edge);
ARRAY_dec(RefEdge);
SET_dec(RefEdge);
ARRAYSET_dec(RefEdge);
ARRAY_dec(ArraysetRefEdge);

#endif
