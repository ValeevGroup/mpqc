
// these provide integrals using the libintv2 routines

#ifndef _chemistry_qc_intv2_obintv2_h
#define _chemistry_qc_intv2_obintv2_h

#include <math/topology/pointbag.h>

#include <chemistry/qc/basis/obint.h>

///////////////////////////////////////////////////////////////////////////

class OneBodyIntv2 : public OneBodyInt {
  private:
    int same_center;

  protected:
    struct struct_centers* c1;
    struct struct_centers* c2;

  public:
    OneBodyIntv2(const RefGaussianBasisSet&, const RefGaussianBasisSet&);
    ~OneBodyIntv2();
};

///////////////////////////////////////////////////////////////////////////

class GaussianOverlapIntv2 : public OneBodyIntv2
{
  public:
    GaussianOverlapIntv2(const RefGaussianBasisSet&,
                         const RefGaussianBasisSet&);
    ~GaussianOverlapIntv2();
    void compute_shell(int,int);
};

class GaussianKineticIntv2 : public OneBodyIntv2
{
  public:
    GaussianKineticIntv2(const RefGaussianBasisSet&,
                         const RefGaussianBasisSet&);
    ~GaussianKineticIntv2();
    void compute_shell(int,int);
};

class GaussianPointChargeIntv2 : public OneBodyIntv2
{
  private:
    int ncharge;
    double** position;
    double* charge;
    RefPointChargeData data_;

  public:
    GaussianPointChargeIntv2(const RefGaussianBasisSet&,
                             const RefGaussianBasisSet&,
                             const RefPointChargeData&);
    ~GaussianPointChargeIntv2();
    void compute_shell(int,int);

    void reinitialize();
};

class GaussianNuclearIntv2 : public GaussianPointChargeIntv2
{
  public:
    GaussianNuclearIntv2(const RefGaussianBasisSet&,
                         const RefGaussianBasisSet&);
    ~GaussianNuclearIntv2();
};

class GaussianEfieldDotVectorIntv2: public OneBodyIntv2
{
  private:
    double *buffer3_; // a larger buffer is needed than that provided
    RefEfieldDotVectorData data_;
  public:
    GaussianEfieldDotVectorIntv2(const RefGaussianBasisSet&,
                                 const RefGaussianBasisSet&,
                                 const RefEfieldDotVectorData&);
    ~GaussianEfieldDotVectorIntv2();
    void compute_shell(int,int);
};

class GaussianDipoleIntv2: public OneBodyIntv2
{
  private:
    RefDipoleData data_;
  public:
    GaussianDipoleIntv2(const RefGaussianBasisSet&,
                        const RefGaussianBasisSet&,
                        const RefDipoleData&);
    ~GaussianDipoleIntv2();
    void compute_shell(int,int);
};

///////////////////////////////////////////////////////////////////////////

class OneBodyDerivIntv2 : public OneBodyDerivInt {
  private:
    int same_center;

    double *intbuf;

  protected:
    struct struct_centers* c1;
    struct struct_centers* c2;

  public:
    OneBodyDerivIntv2(const RefGaussianBasisSet&,
                      const RefGaussianBasisSet&);
    ~OneBodyDerivIntv2();

    void compute_hcore_shell(int center, int ish, int jsh);
    void compute_overlap_shell(int center, int ish, int jsh);
};

#endif
