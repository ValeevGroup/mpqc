
#ifndef _chemistry_qc_intv3_tbintv3_h
#define _chemistry_qc_intv3_tbintv3_h

#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/intv3/int2e.h>

class TwoBodyIntV3 : public TwoBodyInt {
  protected:
    RefInt2eV3 int2ev3_;

  public:
    TwoBodyIntV3(const RefGaussianBasisSet&b1,
                 const RefGaussianBasisSet&b2,
                 const RefGaussianBasisSet&b3,
                 const RefGaussianBasisSet&b4,
                 int storage);
    ~TwoBodyIntV3();

    int log2_shell_bound(int,int,int,int);
    void compute_shell(int,int,int,int);
};

class TwoBodyDerivIntV3 : public TwoBodyDerivInt {
  protected:
    RefInt2eV3 int2ev3_;

  public:
    TwoBodyDerivIntV3(const RefGaussianBasisSet&b1,
                      const RefGaussianBasisSet&b2,
                      const RefGaussianBasisSet&b3,
                      const RefGaussianBasisSet&b4,
                      int storage);
    ~TwoBodyDerivIntV3();

    int log2_shell_bound(int,int,int,int);
    void compute_shell(int,int,int,int,DerivCenters&);
};

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
