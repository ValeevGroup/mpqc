
#ifndef _math_isosurf_tess_h
#define _math_isosurf_tess_h

#include <stdio.h>
#include <math/topology/point.h>

class Tesselation {
  private:
    ArraysetRefPoint& _points;
  protected:
    Tesselation(ArraysetRefPoint&);
  public:
    ~Tesselation();
};

struct dch_tess_struct;
class DirichletTesselation: public Tesselation {
  private:
    dch_tess_struct* _tess;
  public:
    DirichletTesselation(ArraysetRefPoint&);
    ~DirichletTesselation();
    void print();
};

#endif
