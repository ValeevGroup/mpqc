
#ifndef _math_isosurf_vertex_h
#define _math_isosurf_vertex_h

#ifdef __GNUC__
#pragma interface
#endif

extern "C" {
#include <stdlib.h>
#include <math.h>
}
#include <util/ref/ref.h>
#include <util/container/set.h>
#include <math/scmat/matrix.h>
#include <math/isosurf/volume.h>
#include <math/isosurf/implicit.h>

class Vertex: public VRefCount {
  private:
    SCVector3 _point;
    SCVector3 *_normal; // _normal is optional
  public:
    Vertex();
    Vertex(const SCVector3& point,const SCVector3& normal);
    Vertex(const SCVector3& point);
    ~Vertex();
    const SCVector3& point() const { return _point; }
    int has_normal() const { return _normal != 0; }
    const SCVector3& normal() const { return *_normal; }
    void set_point(const SCVector3&p);
    void set_normal(const SCVector3&p);
    operator SCVector3&();

    void print(ostream&o=cout);
};

REF_dec(Vertex);
ARRAY_dec(RefVertex);
SET_dec(RefVertex);
ARRAYSET_dec(RefVertex);
ARRAY_dec(ArraysetRefVertex);

#endif
