
#ifdef __GNUC__
#pragma implementation
#endif

#include <math/topology/pointbag.h>

PointBagElem_double::PointBagElem_double(const Point&p, const double&t):
  obj(t),
  pnt(p)
{
}
PointBagElem_double::PointBagElem_double(const PointBagElem_double&b):
  obj(b.obj),
  pnt(b.pnt)
{
}
PointBagElem_double::~PointBagElem_double()
{
}

Point& PointBagElem_double::point()
{
  return pnt;
}

double& PointBagElem_double::object()
{
  return obj;
}

///////////////////////////////////////////////////////////////////

PointBag_double::PointBag_double()
{
}
PointBag_double::PointBag_double(PointBag_double&b)
{
  for (Pix i=b.first(); i!=0; b.next(i)) {
      add(*((PointBagElem_double*) &b.impl(i)));
    }
}
PointBag_double::~PointBag_double()
{
  for (Pix i=first(); i!=0; next(i)) {
      delete (PointBagElem_double*) impl(i).getptr();
    }
}
