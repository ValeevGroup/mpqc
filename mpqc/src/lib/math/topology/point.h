#ifndef _libQC_point_h
#define _libQC_point_h

// This class implements arbitrary dimensional point using double precision #'s
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include <util/class/class.h>
#include <util/state/state.h>
#include <util/container/ref.h>
#include <util/container/set.h>
#include <util/container/array.h>

class RefSCVector;
class KeyVal;

class Point : virtual public SavableState
{
#   define CLASSNAME Point
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    int dim;
    double *x;
  
  public:
  // Constructors
    Point(int in_dim=3);
    Point(const double *y, int in_dim=3);
    Point(const Point &in_point);
    Point(RefSCVector&);
    Point(double);
    Point(double,double);
    Point(double,double,double);
    Point(StateIn&);
    Point(KeyVal&);
  
  // Destructor
    ~Point(void);
  
    int dimension() const;
    void resize(int dim);

  // Copy to double *
    double *copy(void) const;

    Point& operator=(const Point&);
  
  // Set an element to a particular value
    double &operator[](int i);
    const double &operator[](int i) const;
  
  // Print out a Point
    void print(FILE *fp = stdout);

  // save and restore state
    void save_data_state(StateOut&);
    void restore_data_state(int,StateIn&);

};

DescribedClass_REF_dec(Point);
ARRAY_dec(RefPoint);
SET_dec(RefPoint);
ARRAYSET_dec(RefPoint);

class cart_point {
 private:
  double r[3];
 public:
  cart_point();
  cart_point(const cart_point&p);
  cart_point(const double*);
  ~cart_point();
  double& operator[](int i);
  const double& operator[](int i) const;
  double& x();
  double& y();
  double& z();
  const double& x() const;
  const double& y() const;
  const double& z() const;
};

typedef cart_point Point3;

typedef struct {
  double r;
  double theta;
  double phi;
} sph_point;

#endif


