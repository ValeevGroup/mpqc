
#ifndef _math_isosurf_edge_h
#define _math_isosurf_edge_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/isosurf/vertex.h>

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

#endif
