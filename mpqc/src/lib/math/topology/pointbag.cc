
#ifdef __GNUC__
#pragma implementation
#endif

#include "pointbag.h"

PointBagElem_double::PointBagElem_double(const Point&p, const double&t):
  pnt(p),
  obj(t)
{
}
PointBagElem_double::PointBagElem_double(const PointBagElem_double&b):
  pnt(b.pnt),
  obj(b.obj)
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
      delete (PointBagElem_double*) &impl(i);
    }
}
