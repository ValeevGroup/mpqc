
#ifdef __GNUC__
#pragma implementation
#endif

#include <math/isosurf/vertex.h>


/////////////////////////////////////////////////////////////////////////
// Vertex (a point and a gradient)

REF_def(Vertex);
ARRAY_def(RefVertex);
SET_def(RefVertex);
ARRAYSET_def(RefVertex);
ARRAY_def(ArraysetRefVertex);

Vertex::Vertex(RefSCVector&point,RefSCVector&gradient):
  _point(point),
  _gradient(gradient)
{
}

Vertex::Vertex()
{
}

Vertex::~Vertex()
{
}

RefSCVector
Vertex::gradient()
{
  return _gradient;
}

RefSCVector
Vertex::point()
{
  return _point;
}

void
Vertex::set_point(RefSCVector&p)
{
  _point = p;
}

Vertex::operator RefSCVector()
{
  return _point;
}

RefSCDimension
Vertex::dimension()
{
  return _point.dim();
}

void
Vertex::print(FILE*fp)
{
  int i;
  fprintf(fp, "Vertex:");
  for (i=0; i<_point.dim(); i++)  {
      fprintf(fp," %8.5f", (double)_point[i]);
    }
  for (i=0; i<_gradient.dim(); i++)  {
      fprintf(fp," %8.5f", (double)_gradient[i]);
    }
  fprintf(fp,"\n");
}
