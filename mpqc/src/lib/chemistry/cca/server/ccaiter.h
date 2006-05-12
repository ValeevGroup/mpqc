#ifdef __GNUG__
#pragma interface
#endif

#ifndef _ccaiter_h
#define _ccaiter_h

#include <chemistry/qc/basis/cartiter.h>

namespace MPQC {

class CartesianIterCCA : public sc::CartesianIter {
  int *avec, *bvec, *cvec;
  public:
    CartesianIterCCA(int l) : CartesianIter(l) {}
    void start() {
      bfn_=b_=c_=0;
      a_=l_;
    }
    void next() {
      if (c_ < l_ - a_) {
        b_--;
        c_++;
      }
      else {
        a_--;
        c_ = 0;
        b_ = l_ - a_;
      }
      bfn_++;
    }
    operator int() {
      return (a_ >= 0);
    }
};

}

#endif
