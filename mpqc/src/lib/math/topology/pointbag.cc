
#include "pointbag.h"

PointBagElem_double::PointBagElem_double(Point&p, double&t):
  pnt(p),
  obj(t)
{
}
PointBagElem_double::~PointBagElem_double()
{
}

PointBag_double::PointBag_double()
{
}
PointBag_double::~PointBag_double()
{
  for (Pix i=first(); i!=0; next(i)) {
      delete (PointBagElem_double*) impl(i);
    }
}
