
#ifndef _chemistry_qc_intv3_obintv3_h
#define _chemistry_qc_intv3_obintv3_h

#include <math/topology/pointbag.h>

#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/intv3/int1e.h>

///////////////////////////////////////////////////////////////////////////

class OneBodyIntV3 : public OneBodyInt {
  protected:
    RefInt1eV3 int1ev3_;
    typedef void (Int1eV3::*IntegralFunction)(int,int);
    IntegralFunction intfunc_;
  public:
    OneBodyIntV3(const RefGaussianBasisSet&, const RefGaussianBasisSet&,
                 IntegralFunction);
    ~OneBodyIntV3();
    void compute_shell(int,int);
};

class PointChargeIntV3 : public OneBodyInt
{
  protected:
    RefInt1eV3 int1ev3_;
    int ncharge;
    double** position;
    double* charge;
    RefPointChargeData data_;
  public:
    PointChargeIntV3(const RefGaussianBasisSet&,
                     const RefGaussianBasisSet&,
                     const RefPointChargeData&);
    ~PointChargeIntV3();
    void compute_shell(int,int);

    void reinitialize();
};

class EfieldDotVectorIntV3: public OneBodyInt
{
  protected:
    RefInt1eV3 int1ev3_;
    RefEfieldDotVectorData data_;
  public:
    EfieldDotVectorIntV3(const RefGaussianBasisSet&,
                         const RefGaussianBasisSet&,
                         const RefEfieldDotVectorData&);
    ~EfieldDotVectorIntV3();
    void compute_shell(int,int);
};

class DipoleIntV3: public OneBodyInt
{
  protected:
    RefInt1eV3 int1ev3_;
    RefDipoleData data_;
  public:
    DipoleIntV3(const RefGaussianBasisSet&,
                const RefGaussianBasisSet&,
                const RefDipoleData&);
    ~DipoleIntV3();
    void compute_shell(int,int);
};

///////////////////////////////////////////////////////////////////////////

class OneBodyDerivIntV3 : public OneBodyDerivInt {
  protected:
    RefInt1eV3 int1ev3_;
    typedef void (Int1eV3::*IntegralFunction)(int,int,int,int);
    IntegralFunction intfunc_;
  public:
    OneBodyDerivIntV3(const RefGaussianBasisSet&,
                      const RefGaussianBasisSet&,
                      IntegralFunction);
    ~OneBodyDerivIntV3();
    void compute_shell(int,int,DerivCenters&);
};

#endif
