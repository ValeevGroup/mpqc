
#ifndef _math_isosurf_tricoef_h
#define _math_isosurf_tricoef_h

#include <util/container/ref.h>

class TriInterpCoefKey {
  private:
    int order_;
    double L2_;
    double L3_;
  public:
    TriInterpCoefKey(int order, double L2, double L3):
      order_(order), L2_(L2), L3_(L3) {}
    int order() const { return order_; }
    double L1() const { return 1.0 - L2_ - L3_; }
    double L2() const { return L2_; }
    double L3() const { return L3_; }
    int cmp(const TriInterpCoefKey&t) const {
        if (order_ < t.order_) return -1;
        if (order_ > t.order_) return 1;
        if (L2_ < t.L2_) return -1;
        if (L2_ > t.L2_) return 1;
        if (L3_ < t.L3_) return -1;
        if (L3_ > t.L3_) return 1;
        return 0;
      }
};

#define TriInterpCoefKeyEQ(k1,k2) ((k1).cmp(k2)==0)
#define TriInterpCoefKeyCMP(k1,k2) ((k1).cmp(k2))

class TriInterpCoef: public VRefCount {
    double *coef_;
    double *rderiv_;
    double *sderiv_;
  public:
    TriInterpCoef(const TriInterpCoefKey& key);
    ~TriInterpCoef();
    double& coef(int i, int j, int k) {return coef_[ijk_to_index(i,j,k)];}
    double& rderiv(int i, int j, int k) {return rderiv_[ijk_to_index(i,j,k)];}
    double& sderiv(int i, int j, int k) {return sderiv_[ijk_to_index(i,j,k)];}

    static int
    ijk_to_index(int i, int j, int k)
    {
      int n = i + j + k;
      int ir = n - i;
      return (ir*(ir+1)>>1) + j;
    }

    static int
    order_to_nvertex(int order)
    {
      return ((order+1)*(order+2)>>1);
    }
};

Ref_declare(TriInterpCoef);

#endif
