
// these provide integrals using Justin Fermann's cints routines

#ifndef _chemistry_qc_cints_integraljf_h
#define _chemistry_qc_cints_integraljf_h

#include <chemistry/qc/integral/integral.h>
#include <math/topology/pointbag.h>

class OneBodyIntJF : public OneBodyInt
{
  protected:
    static int maxfact;
    
    int maxam;
    int *lci;
    unsigned int *num_ser;
    
    double *df;
    
    void init();
    double f_n(int, int, int, double, double);
    double overlap_int(double, int, int, int, double,
                       double, int, int, int, double,
                       double, double, Point&, Point&, int);

  public:
    OneBodyIntJF(const RefGaussianBasisSet&, OneBodyIntIter* =0);
    OneBodyIntJF(const RefGaussianBasisSet&, const RefGaussianBasisSet&,
                 OneBodyIntIter* =0);
    ~OneBodyIntJF();
};

class GaussianOverlapIntJF : public OneBodyIntJF
{
  public:
    GaussianOverlapIntJF(const RefGaussianBasisSet&, OneBodyIntIter* =0);
    GaussianOverlapIntJF(const RefGaussianBasisSet&,
                         const RefGaussianBasisSet&,
                         OneBodyIntIter* =0);
    ~GaussianOverlapIntJF();
    void compute_shell(int,int,double*);
};

class GaussianKineticIntJF : public OneBodyIntJF
{
  protected:
    double ke_int(double, int, int, int, double,
                  double, int, int, int, double,
                  double, double, Point&, Point&, int);
  public:
    GaussianKineticIntJF(const RefGaussianBasisSet&, OneBodyIntIter* =0);
    GaussianKineticIntJF(const RefGaussianBasisSet&,
                         const RefGaussianBasisSet&,
                         OneBodyIntIter* =0);
    ~GaussianKineticIntJF();
    void compute_shell(int,int,double*);
};

#if 0
class GaussianPointChargeIntJF : public OneBodyIntJF
{
  private:
    int ncharge;
    double** position;
    double* charge;

    void init(PointBag_double*);
    
  public:
    GaussianPointChargeIntv2(PointBag_double*, const RefGaussianBasisSet&,
                             OneBodyIntIter* =0);
    GaussianPointChargeIntv2(PointBag_double*, const RefGaussianBasisSet&,
                             const RefGaussianBasisSet&,
                             OneBodyIntIter* =0);
    ~GaussianPointChargeIntv2();
    void compute_shell(int,int,double*);
};

class GaussianNuclearIntJF : public GaussianPointChargeIntJF
{
  public:
    GaussianNuclearIntv2(const RefGaussianBasisSet&, OneBodyIntIter* =0);
    GaussianNuclearIntv2(PointBag_double *charges,
                         const RefGaussianBasisSet&,
                         const RefGaussianBasisSet&,
                         OneBodyIntIter* =0);

    ~GaussianNuclearIntv2();
};
#endif

#endif
