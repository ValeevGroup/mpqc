
#ifdef __GNUC__
#pragma interface
#endif

#ifdef INLINE_FUNCTIONS
#define INLINE inline
#else
#define INLINE
#endif

INLINE void
SCVector3::spherical_coord(double theta, double phi,
                           double r)
{
    _v[0]=r*sin(theta)*cos(phi);
    _v[1]=r*sin(theta)*sin(phi);
    _v[2]=r*cos(theta);

}

INLINE  double
SCVector3::dist(const SCVector3 &s) const
{
    return sqrt((_v[0]-s._v[0])*(_v[0]-s._v[0])+
                (_v[1]-s._v[1])*(_v[1]-s._v[1])+
                (_v[2]-s._v[2])*(_v[2]-s._v[2]));
}

#undef INLINE
