
#ifndef _math_isosurf_volume_h
#define _math_isosurf_volume_h

#include <math/opt/nlp.h>
#include <util/container/ref.h>
#include <math/nihmatrix/nihmatrix.h>
#include <math/topology/point.h>

class Volume: public NLP2 {
#   define CLASSNAME Volume
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  private:
    double _interp_acc;

    Resultdouble& _value;
    ResultColumnVector& _gradient;
    ResultSymmetricMatrix& _hessian;
  protected:
    void set_value(double);
    void set_gradient(ColumnVector&);
    void set_gradient(DVector&);
    void set_hessian(SymmetricMatrix&);

    inline double& interpolation_accuracy() { return _interp_acc; };

    virtual void compute() = 0;

    virtual void failure(const char*);
  public:
    Volume(int dimension);
    ~Volume();
    void SetX(const Point& x);
    inline void SetX(RefPoint& x) { SetX(*(x.pointer())); }
    inline void SetX(ColumnVector& x) { NLP0::SetX(x); }
    void GetX(const Point& x);
    inline const ColumnVector& GetX() { return NLP0::GetX(); }

    void gradient(DVector&); // calls gradient() and converts it to DVector
    inline void gradient(RefPoint&p, DVector&g)
    {
      SetX(p); gradient(g);
    }

    // find the corners of a bounding box which approximately
    // contains all points with a value between valuemin and valuemax
    // the result must satisfy p1[i] < p2[i]
    virtual void boundingbox(double valuemin,
                             double valuemax,
                             Point& p1, Point& p2) = 0;
    virtual void pointset(double resolution,
                          double valuemin,
                          double valuemax,
                          SetRefPoint& points);

    virtual RefPoint interpolate(RefPoint,RefPoint,double value);
    virtual RefPoint solve(RefPoint,DVector&,double value);
    
    int do_value(int);
    int do_gradient(int);
    int do_hessian(int);
    int do_value();
    int do_gradient();
    int do_hessian();

    virtual void Eval();
    virtual double EvalF();
    virtual ColumnVector EvalG();
    virtual SymmetricMatrix EvalH();

    virtual double value();
    virtual inline double value(RefPoint&p) { SetX(p); return value(); };
    virtual const ColumnVector& gradient();
    virtual const SymmetricMatrix& hessian();
};

SavableState_REF_dec(Volume);
ARRAY_dec(RefVolume);
SET_dec(RefVolume);
ARRAYSET_dec(RefVolume);

#endif
