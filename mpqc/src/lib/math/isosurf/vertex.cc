
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

Vertex::Vertex(const SCVector3&point):
  _point(point),
  _normal(0)
{
}

Vertex::Vertex(const SCVector3&point,const SCVector3&normal):
  _point(point),
  _normal(new SCVector3(normal))
{
}

Vertex::Vertex():
  _normal(0)
{
}

Vertex::~Vertex()
{
  if (_normal) delete _normal;
}

void
Vertex::set_point(const SCVector3&p)
{
  _point = p;
}

void
Vertex::set_normal(const SCVector3&p)
{
  if (_normal) {
      *_normal = p;
    }
  else {
      _normal = new SCVector3(p);
    }
}

Vertex::operator SCVector3&()
{
  return _point;
}

void
Vertex::print(FILE*fp)
{
  int i;
  fprintf(fp, "Vertex:");
  for (i=0; i<3; i++)  {
      fprintf(fp," %8.5f", _point[i]);
    }
  if (_normal) {
      for (i=0; i<3; i++)  {
          fprintf(fp," %8.5f", normal()[i]);
        }
    }
  fprintf(fp,"\n");
}
