
#ifndef _math_isosurf_isosurf_h
#define _math_isosurf_isosurf_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/isosurf/surf.h>

class IsosurfaceGen {
  protected:
    double _resolution;
  public:
    IsosurfaceGen();
    virtual ~IsosurfaceGen();
    virtual void isosurface(double value,
                            TriangulatedSurface& surf) = 0;
    virtual void set_resolution(double);
};

class ImplicitSurfacePolygonizer: public IsosurfaceGen {
  private:
    // These static data and members are used to interface to the
    // implicit.c routine provided in Graphics Gems IV.
    static ImplicitSurfacePolygonizer* current;
    static int add_triangle_to_current(int,int,int,VERTICES);
    static double value_of_current(double x, double y, double z);
  protected:
    RefVolume _volume;

    ArraysetRefVertex  _tmp_vertices;
    TriangulatedSurface* _surf;
    double _value;
  public:
    ImplicitSurfacePolygonizer(const RefVolume&);
    virtual ~ImplicitSurfacePolygonizer();
    virtual void isosurface(double value,
                            TriangulatedSurface& surf);
};  

#endif
