
#ifdef INLINE_FUNCTIONS
#define INLINE inline
#else
#define INLINE
#endif

INLINE int PointBag_double::length()
{
  return impl.length();
}

INLINE void PointBag_double::add(const Point& p, double x)
{
  impl.add((void*)new PointBagElem_double(p,x));
}

INLINE void PointBag_double::add(const PointBagElem_double& e)
{
  impl.add((void*)new PointBagElem_double(e));
}

INLINE int PointBag_double::owns(Pix i)
{
  return impl.owns(i);
}

INLINE double& PointBag_double::operator()(Pix i)
{
  return ((PointBagElem_double*)&impl(i))->object();
}

INLINE double& PointBag_double::get(Pix i)
{
  return operator()(i);
}

INLINE Point& PointBag_double::point(Pix i)
{
  return ((PointBagElem_double*)&impl(i))->point();
}

INLINE Pix PointBag_double::first()
{
  return impl.first();
}

INLINE void PointBag_double::next(Pix&i)
{
  impl.next(i);
}

#undef INLINE
