
// these provide integrals using Justin Fermann's cints routines

#ifndef _chemistry_qc_cints_integraljf_h
#define _chemistry_qc_cints_integraljf_h

#include <chemistry/qc/basis/integral.h>
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
    double overlap_intx(int, int, double, double, double);
    double overlap_int(int, int, int, int, int, int,
                       double[3], double[3], double);

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
    double ke_int(int, int, int, int, int, int, double[3], double[3],
                  double, double, double);

  public:
    GaussianKineticIntJF(const RefGaussianBasisSet&, OneBodyIntIter* =0);
    GaussianKineticIntJF(const RefGaussianBasisSet&,
                         const RefGaussianBasisSet&,
                         OneBodyIntIter* =0);
    ~GaussianKineticIntJF();
    void compute_shell(int,int,double*);
};

class GaussianPointChargeIntJF : public OneBodyIntJF
{
  private:
    int ncharge;
    double *charge;
    double **position;

    double ****vint;
    char ***vint_done;

    void init(PointBag_double*);
    void init_peint(int, int, int, int, double, double, double,
                    double[3], double[3], double[3]);
    double * do_pe_doublet(double[3][3], double[], int, int, int, int, double);
    void calc_f(double[], int, double);
    
  public:
    GaussianPointChargeIntJF(PointBag_double*, const RefGaussianBasisSet&,
                             OneBodyIntIter* =0);
    GaussianPointChargeIntJF(PointBag_double*, const RefGaussianBasisSet&,
                             const RefGaussianBasisSet&,
                             OneBodyIntIter* =0);
    ~GaussianPointChargeIntJF();
    void compute_shell(int,int,double*);
};

class GaussianNuclearIntJF : public GaussianPointChargeIntJF
{
  public:
    GaussianNuclearIntJF(const RefGaussianBasisSet&, OneBodyIntIter* =0);
    GaussianNuclearIntJF(PointBag_double *charges,
                         const RefGaussianBasisSet&,
                         const RefGaussianBasisSet&,
                         OneBodyIntIter* =0);

    ~GaussianNuclearIntJF();
};

#endif
