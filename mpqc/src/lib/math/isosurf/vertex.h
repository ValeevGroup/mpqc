
#ifndef _math_isosurf_vertex_h
#define _math_isosurf_vertex_h

#ifdef __GNUC__
#pragma interface
#endif

extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
}
#include <util/container/ref.h>
#include <util/container/set.h>
#include <math/scmat/matrix.h>
#include <math/isosurf/volume.h>
#include <math/isosurf/implicit.h>

class Vertex: public RefCount {
  private:
    RefSCVector _point;
    RefSCVector _gradient;
  public:
    Vertex();
    Vertex(RefSCVector&point,RefSCVector&gradient);
    ~Vertex();
    RefSCVector gradient();
    RefSCVector point();
    void set_point(RefSCVector&p);
    operator RefSCVector();
    RefSCDimension dimension();

    void print(FILE*fp = stdout);
};

REF_dec(Vertex);
ARRAY_dec(RefVertex);
SET_dec(RefVertex);
ARRAYSET_dec(RefVertex);
ARRAY_dec(ArraysetRefVertex);

#endif
