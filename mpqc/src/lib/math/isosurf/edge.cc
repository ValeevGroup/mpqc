
#ifdef __GNUC__
#pragma implementation
#endif

#include <math/isosurf/edge.h>

/////////////////////////////////////////////////////////////////////////
// Edge

REF_def(Edge);
ARRAY_def(RefEdge);
SET_def(RefEdge);
ARRAYSET_def(RefEdge);
ARRAY_def(ArraysetRefEdge);

Edge::~Edge()
{
}

double Edge::length()
{
  RefSCVector A(_vertices[0]->point());
  RefSCVector B(_vertices[1]->point());
  RefSCVector BA = B - A;
  return sqrt(BA.dot(BA));
}

void Edge::add_vertices(SetRefVertex&set)
{
  set.add(_vertices[0]);
  set.add(_vertices[1]);
}

/////////////////////////////////////////////////////////////////////////
// Edge4

REF_def(Edge4);
ARRAY_def(RefEdge4);
SET_def(RefEdge4);
ARRAYSET_def(RefEdge4);
ARRAY_def(ArraysetRefEdge4);

Edge4::~Edge4()
{
}

Edge4::Edge4(RefVertex p1,RefVertex p2,const RefVolume&vol,double isovalue):
  Edge(p1,p2)
{
  //// find the initial guess for the interior vertices
  // Note: interior vertex 0 is next to vertex 0
  RefVertex p[2];
  p[0] = p1;
  p[1] = p2;
  RefSCVector pv[2];
  RefSCVector grad[2];

  int i;
  for (i=0; i<2; i++) {
      pv[i] = p[i]->point();
      grad[i] = p[i]->gradient();
    }

  for (i=0; i<2; i++) {
      RefSCVector interpv;
      interpv = ((2.0*pv[i])+pv[(i==1)?0:1])*(1.0/3.0);
      RefSCVector start(interpv.dim());
      start.assign(interpv);
      RefSCVector interpgrad;
      interpgrad = ((2.0*grad[i])+grad[(i==1)?0:1])*(1.0/3.0);
      RefSCVector newpoint = vol->solve(start,interpgrad,isovalue);
      vol->set_x(newpoint);
      interpgrad = vol->gradient().copy();
      _interiorvertices[i] = new Vertex(newpoint,interpgrad);
    }

}

// Returns a length assuming linearity for now.
double
Edge4::length()
{
  //fprintf(stderr,"Edge4::length(): not implemented\n");
  //abort();
  return Edge::length();
}
