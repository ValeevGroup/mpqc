
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

Edge::Edge(const RefVertex &p1,
     const RefVertex &p2,
     const RefVertex &p3)
{
  _order = 2;
  _vertices = new RefVertex[3];
  _vertices[0]=p1; _vertices[1]=p2; _vertices[2]=p3;
}

Edge::Edge(const RefVertex &p1,
     const RefVertex &p2,
     const RefVertex &p3,
     const RefVertex &p4)
{
  _order = 3;
  _vertices = new RefVertex[4];
  _vertices[0]=p1; _vertices[1]=p2; _vertices[2]=p3; _vertices[3]=p4;
}

Edge::~Edge()
{
  delete[] _vertices;
}

double Edge::straight_length()
{
  SCVector3 BA = vertex(1)->point() - vertex(0)->point();
  return BA.norm();
}

void Edge::add_vertices(SetRefVertex&set)
{
  set.add(_vertices[0]);
  set.add(_vertices[_order]);
}

void
Edge::set_order(int order, const RefVolume&vol,double isovalue)
{
  RefVertex *newvertices = new RefVertex[order+1];
  newvertices[0] = vertex(0);
  newvertices[order] = vertex(1);
  delete[] _vertices;
  _vertices = newvertices;
  _order = order;

  SCVector3 pv[2];
  SCVector3 norm[2];

  int i;
  for (i=0; i<2; i++) {
      pv[i] = vertex(i)->point();
      norm[i] = vertex(i)->normal();
    }

  for (i=1; i<_order; i++) {
      double I = (1.0*i)/order;
      double J = (1.0*(_order - i))/order;
      SCVector3 interpv = I * pv[0] + J * pv[1];
      SCVector3 start(interpv);
      SCVector3 interpnorm = I * norm[0] + J * norm[1];
      SCVector3 newpoint;
      vol->solve(start,interpnorm,isovalue,newpoint);
      vol->set_x(newpoint);
      if (vol->gradient_implemented()) {
          vol->get_gradient(interpnorm);
          interpnorm.normalize();
        }
      _vertices[i] = new Vertex(newpoint,interpnorm);
    }

}
