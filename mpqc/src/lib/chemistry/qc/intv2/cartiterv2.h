
// these provide integrals using the libintv2 routines

#ifndef _chemistry_qc_intv2_cartiterv2_h
#define _chemistry_qc_intv2_cartiterv2_h

#include <chemistry/qc/basis/cartiter.h>

class CartesianIterV2 : public CartesianIter {
  public:
    CartesianIterV2(int l) : CartesianIter(l) {}

    void start() {
      bfn_=a_=c_=0;
      b_=l_;
    }

    void next() {
      if (c_<l_-a_)
        c_++;
      else {
        c_=0;
        a_++;
      }
      bfn_++;
      b_ = l_-a_-c_;
    }
    
    operator int() {
      return (a_ <= l_);
    }
};

class RedundantCartesianIterV2 : public RedundantCartesianIter {
  public:
    RedundantCartesianIterV2(int l) : RedundantCartesianIter(l) {}

    int bfn() {
      int i = a();
      int j = b();
      int am = l();
      return (((((((am)+1)<<1)-(i))*((i)+1))>>1)-(j)-1);
    }
};

class RedundantCartesianSubIterV2 : public RedundantCartesianSubIter {
  public:
    RedundantCartesianSubIterV2(int l) : RedundantCartesianSubIter(l) {}

    int bfn() {
      int i = a();
      int j = b();
      int am = l();
      return (((((((am)+1)<<1)-(i))*((i)+1))>>1)-(j)-1);
    }
};

#endif
