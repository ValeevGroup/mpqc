
#if defined(__GNUC__)
#pragma interface
#endif

#ifndef _math_symmetry_tform_h
#define _math_symmetry_tform_h

#include <math/symmetry/pointgrp.h>

class SymSphericalTransform {
  public:
    class Component {
      private:
        double coef_;
        int a_, b_, c_, cartindex_, pureindex_;
      public:
        int a() const { return a_; }
        int b() const { return b_; }
        int c() const { return c_; }
        int cartindex() const { return cartindex_; }
        int pureindex() const { return pureindex_; }
        double coef() const { return coef_; }
        void init(int a, int b, int c, double coef, int pureindex);
    };
  protected:
    int n_;
    int l_;
    Component *components_;
    SymSphericalTransform();
  public:
    SymSphericalTransform(int l);
    ~SymSphericalTransform();
    int cartindex(int i) const { return components_[i].cartindex(); }
    int pureindex(int i) const { return components_[i].pureindex(); }
    double coef(int i) const { return components_[i].coef(); }
    int a(int i) const { return components_[i].a(); }
    int b(int i) const { return components_[i].b(); }
    int c(int i) const { return components_[i].c(); }
    int l() const { return l_; }
    int n() const { return n_; }

    void add(int a, int b, int c, double coef, int pureindex);
};

// The inverse transforms
class ISymSphericalTransform: public SymSphericalTransform {
  public:
    ISymSphericalTransform(int l);
};

class SymSphericalTransformIter {
  private:
    const SymSphericalTransform *transform_;

    int i_;
  public:
    SymSphericalTransformIter(SymSphericalTransform*);
    SymSphericalTransformIter(int l, int inverse = 0);
    void begin() { i_ = 0; }
    void start() { begin(); }
    void next() { i_++; }
    int ready() { return i_ < transform_->n(); }
    operator int() { return ready(); }
    int l() { return transform_->l(); }
    int cartindex() { return transform_->cartindex(i_); }
    int pureindex() { return transform_->pureindex(i_); }
    int bfn() { return pureindex(); }
    double coef() { return transform_->coef(i_); }
    int a() { return transform_->a(i_); }
    int b() { return transform_->b(i_); }
    int c() { return transform_->c(i_); }
    int l(int i) { return i?(i==1?b():c()):a(); }
    int n() { return 2*l() + 1; }
};

class SymCartesianIter
{
  private:
    int a_;
    int b_;
    int c_;
    int l_;
    int bfn_;
  public:
    SymCartesianIter(int l);
    ~SymCartesianIter();
    void start();
    void next();
    operator int();
    int a();
    int b();
    int c();
    int l(int i);
    int l();
    int bfn();
    int n();
};

class SymRedundantCartesianIter {
  private:
    int done_;
    int l_;
    int *axis_;
  public:
    SymRedundantCartesianIter(int l);
    ~SymRedundantCartesianIter();
    void start();
    void next();
    operator int();
    int bfn();
    int l(int i);
    int l() { return l_; }
    int a() { return l(0); }
    int b() { return l(1); }
    int c() { return l(2); }
    int axis(int i) { return axis_[i]; }
};

// Like RedundantCartesianIter, except a, b, and c are fixed to a given value
class SymRedundantCartesianSubIter {
  private:
    int done_;
    int l_;
    int e_[3];
    int *axis_;

    void advance();
    int valid();
  public:
    SymRedundantCartesianSubIter(int l);
    ~SymRedundantCartesianSubIter();
    void start(int a, int b, int c);
    void next();
    operator int();
    int bfn();
    int l(int i) { return e_[i]; }
    int l() { return l_; }
    int a() { return e_[0]; }
    int b() { return e_[1]; }
    int c() { return e_[2]; }
    int axis(int i) { return axis_[i]; }
};

class SymRotation {
  private:
    int _n;
    int _am;
    double **r;
    
    void done();

  public:
    void init(int a, SymmetryOperation&so);
    void init_pure(int a, SymmetryOperation&so);
    
    SymRotation(int a, SymmetryOperation& so, int pure = 0);
    ~SymRotation();

    int am() const { return _am; }
    int dim() const { return _n; }
    
    double& operator()(int i, int j) { return r[i][j]; }
    double* operator[](int i) { return r[i]; }
    
    void print() const;
    
    double trace() const;
};

#endif
