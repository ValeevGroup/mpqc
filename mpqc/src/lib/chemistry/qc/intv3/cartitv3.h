
#ifndef _chemistry_qc_intv3_cartitv3_h
#define _chemistry_qc_intv3_cartitv3_h

#include <chemistry/qc/basis/cartiter.h>

class CartesianIterV3 : public CartesianIter {
  public:
    CartesianIterV3(int l) : CartesianIter(l) {}

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

class RedundantCartesianIterV3 : public RedundantCartesianIter {
  public:
    RedundantCartesianIterV3(int l) : RedundantCartesianIter(l) {}

    int bfn() {
      int i = a();
      int j = b();
      int am = l();
      return (((((((am)+1)<<1)-(i))*((i)+1))>>1)-(j)-1);
    }
};

class RedundantCartesianSubIterV3 : public RedundantCartesianSubIter {
  public:
    RedundantCartesianSubIterV3(int l) : RedundantCartesianSubIter(l) {}

    int bfn() {
      int i = a();
      int j = b();
      int am = l();
      return (((((((am)+1)<<1)-(i))*((i)+1))>>1)-(j)-1);
    }
};

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
