
#ifndef _chemistry_qc_basis_cartiter_h
#define _chemistry_qc_basis_cartiter_h

#ifdef __GNUC__
#pragma interface
#endif

class CartesianIter {
  protected:
    int a_;
    int b_;
    int c_;
    int l_;
    int bfn_;

  public:
    CartesianIter(int l);
    virtual ~CartesianIter();

    virtual void start() =0;
    virtual void next() =0;
    virtual operator int() =0;

    int n();
    int a();
    int b();
    int c();
    int l();
    int l(int i);
    int bfn();
};

class RedundantCartesianIter {
  private:
    int done_;
    int l_;
    int *axis_;

  public:
    RedundantCartesianIter(int l);
    virtual ~RedundantCartesianIter();

    virtual int bfn() =0;

    void start();
    void next();
    operator int();

    int a();
    int b();
    int c();
    int l();
    int l(int i);
    int axis(int i);
};

// Like RedundantCartesianIter, except a, b, and c are fixed to a given value
class RedundantCartesianSubIter {
  private:
    int done_;
    int l_;
    int e_[3];
    int *axis_;

    void advance();
    int valid();

  public:
    RedundantCartesianSubIter(int l);
    ~RedundantCartesianSubIter();

    virtual int bfn() =0;

    void start(int a, int b, int c);
    void next();
    operator int();

    int l();
    int a();
    int b();
    int c();
    int l(int i);
    int axis(int i);
};

#ifdef INLINE_FUNCTIONS
#include <chemistry/qc/basis/cartiter_i.h>
#endif
  
#endif
