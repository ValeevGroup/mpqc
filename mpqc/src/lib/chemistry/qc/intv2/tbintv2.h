
#ifndef _chemistry_qc_intv2_tbintv2_h
#define _chemistry_qc_intv2_tbintv2_h

#include <chemistry/qc/basis/tbint.h>

class TwoBodyIntV2 : public TwoBodyInt {
  private:
    int same_center;
    struct struct_centers* c1;
    struct struct_centers* c2;
    struct struct_centers* c3;
    struct struct_centers* c4;

    double *intbuf;
    
    void init();
    
  public:
    TwoBodyIntV2(const RefGaussianBasisSet&b);
    TwoBodyIntV2(const RefGaussianBasisSet&b1,
                 const RefGaussianBasisSet&b2,
                 const RefGaussianBasisSet&b3,
                 const RefGaussianBasisSet&b4);
    virtual ~TwoBodyIntV2();

    void compute_shell(int,int,int,int);
};

class TwoBodyDerivIntV2 : public TwoBodyDerivInt {
  private:
    int same_center;
    struct struct_centers* c1;
    struct struct_centers* c2;
    struct struct_centers* c3;
    struct struct_centers* c4;

    double *intbuf;
    
    void init();
    
  public:
    TwoBodyDerivIntV2(const RefGaussianBasisSet&b);
    TwoBodyDerivIntV2(const RefGaussianBasisSet&b1,
                      const RefGaussianBasisSet&b2,
                      const RefGaussianBasisSet&b3,
                      const RefGaussianBasisSet&b4);
    virtual ~TwoBodyDerivIntV2();

    void compute_shell(int,int,int,int);
};

#endif
