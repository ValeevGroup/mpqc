
// these provide integrals using the libintv2 routines

#ifndef _chemistry_qc_integral_integralv2_h
#define _chemistry_qc_integral_integralv2_h

#include <chemistry/qc/integral/integral.h>
#include <math/topology/pointbag.h>

class OneBodyIntv2:
  public OneBodyInt
{
 protected:
  struct struct_centers* c;
 public:
  OneBodyIntv2(const GaussianBasisSet*,const Molecule*);
  ~OneBodyIntv2();
};

class GaussianOverlapIntv2:
  public OneBodyIntv2
{
 private:
 public:
  GaussianOverlapIntv2(const GaussianBasisSet*);
  ~GaussianOverlapIntv2();
  void compute_shell(int,int,double*);
};

class GaussianKineticIntv2:
  public OneBodyIntv2
{
 private:
 public:
  GaussianKineticIntv2(const GaussianBasisSet*,const Molecule*);
  ~GaussianKineticIntv2();
  void compute_shell(int,int,double*);
};

class GaussianPointChargeIntv2:
  public OneBodyIntv2
{
 private:
  int ncharge;
  double** position;
  double* charge;
 public:
  GaussianPointChargeIntv2(PointBag_double*,const GaussianBasisSet*,const Molecule*);
  ~GaussianPointChargeIntv2();
  void compute_shell(int,int,double*);
};

class GaussianNuclearIntv2:
  public GaussianPointChargeIntv2
{
 public:
  GaussianNuclearIntv2(const GaussianBasisSet*,const Molecule*);
  ~GaussianNuclearIntv2();
};

#endif
