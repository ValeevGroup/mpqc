
#ifdef __GNUC__
#pragma interface
#endif

#ifdef INLINE_FUNCTIONS
#define INLINE inline
#else
#define INLINE
#endif

INLINE
SCVector3::SCVector3(double x[3])
{
  _v[0] = x[0];
  _v[1] = x[1];
  _v[2] = x[2];
};

INLINE
SCVector3::SCVector3(double x,double y,double z)
{
  _v[0] = x;
  _v[1] = y;
  _v[2] = z;
};

INLINE
SCVector3::SCVector3(const SCVector3&p)
{
  _v[0] = p[0];
  _v[1] = p[1];
  _v[2] = p[2];
};

INLINE double
SCVector3::norm() const
{
  return sqrt(this->dot(*this));
};

INLINE const double*
SCVector3::data() const
{             
  return _v;
}             

INLINE double&
SCVector3::x()
{             
  return _v[0];
}             
              
INLINE double&
SCVector3::y()
{             
  return _v[1];
}             
              
INLINE double&
SCVector3::z()
{
  return _v[2];
}

INLINE const double&
SCVector3::x() const
{
  return _v[0];
}

INLINE const double&
SCVector3::y() const
{
  return _v[1];
}

INLINE const double&
SCVector3::z() const
{
  return _v[2];
}

#undef INLINE
