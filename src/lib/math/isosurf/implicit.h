
#ifndef _math_isosurf_implicit_h
#define _math_isosurf_implicit_h

namespace sc {
  namespace detail {

    typedef struct point { /* a three-dimensional point */
        double x, y, z; /* its coordinates */
    } POINT;

    typedef struct vertex { /* surface vertex */
        POINT position, normal; /* position and surface normal */
    } VERTEX;

    typedef struct vertices { /* list of vertices in polygonization */
        int count, max; /* # vertices, max # allowed */
        VERTEX *ptr; /* dynamically allocated */
    } VERTICES;

  }
}

#define TET	0  /* use tetrahedral decomposition */
#define NOTET	1  /* no tetrahedral decomposition  */

extern "C" {
    char * polygonize(double(*function)(double,double,double),
                      double size, int bounds,
                      double x, double y, double z,
                      int(*triproc)(int,int,int,sc::detail::VERTICES), int mode);
}

#endif
