
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

Vertex::Vertex(RefSCVector&point,RefSCVector&normal):
  _point(point),
  _normal(normal)
{
}

Vertex::Vertex()
{
}

Vertex::~Vertex()
{
}

void
Vertex::set_point(const RefSCVector&p)
{
  _point = p;
}

void
Vertex::set_normal(const RefSCVector&p)
{
  _normal = p;
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
  if (_normal.nonnull()) {
      for (i=0; i<_normal.dim(); i++)  {
          fprintf(fp," %8.5f", (double)_normal[i]);
        }
    }
  fprintf(fp,"\n");
}
