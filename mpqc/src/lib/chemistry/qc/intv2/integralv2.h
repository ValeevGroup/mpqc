
// these provide integrals using the libintv2 routines

#ifndef _chemistry_qc_intv2_integralv2_h
#define _chemistry_qc_intv2_integralv2_h

#include <math/topology/pointbag.h>

#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/cartiter.h>
#include <chemistry/qc/basis/transform.h>

#include <chemistry/qc/intv2/int_macros.h>
#include <chemistry/qc/intv2/atoms.h>

///////////////////////////////////////////////////////////////////////////

class CartesianIterV2 : public CartesianIter {
  public:
    CartesianIterV2(int l) : CartesianIter(l) {}

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

class RedundantCartesianIterV2 : public RedundantCartesianIter {
  public:
    RedundantCartesianIterV2(int l) : RedundantCartesianIter(l) {}

    int bfn() {
      int i = a();
      int j = b();
      int am = l();
      return (((((((am)+1)<<1)-(i))*((i)+1))>>1)-(j)-1);
    }
};

class RedundantCartesianSubIterV2 : public RedundantCartesianSubIter {
  public:
    RedundantCartesianSubIterV2(int l) : RedundantCartesianSubIter(l) {}

    int bfn() {
      int i = a();
      int j = b();
      int am = l();
      return (((((((am)+1)<<1)-(i))*((i)+1))>>1)-(j)-1);
    }
};

///////////////////////////////////////////////////////////////////////////

class SphericalTransformComponentV2 : public SphericalTransformComponent {
  public:
    void init(int a, int b, int c, double coef, int pureindex) {
      a_ = a;
      b_ = b;
      c_ = c;
      coef_ = coef;
      pureindex_ = pureindex;
      cartindex_ = INT_CARTINDEX(a+b+c,a,b);
    }
};

class SphericalTransformV2 : public SphericalTransform {
  public:
    SphericalTransformV2(int l) {
      n_=0;
      l_=l;
      components_=0;
      init();
    }

    SphericalTransformComponent * new_components() {
      return new SphericalTransformComponentV2[n_+1];
    }
};

class ISphericalTransformV2 : public ISphericalTransform {
  public:
    ISphericalTransformV2(int l) {
      n_ = 0;
      l_ = l;
      components_ = 0;
      init();
    }

    SphericalTransformComponent * new_components() {
      return new SphericalTransformComponentV2[n_+1];
    }
};

class SphericalTransformIterV2 : public SphericalTransformIter {
  public:
    SphericalTransformIterV2(int l, int inverse=0);
};

///////////////////////////////////////////////////////////////////////////

class OneBodyIntv2 : public OneBodyInt
{
  private:
    int same_center;
  protected:
    struct struct_centers* c1;
    struct struct_centers* c2;
  public:
    OneBodyIntv2(const RefGaussianBasisSet&, OneBodyIntIter* =0);
    OneBodyIntv2(const RefGaussianBasisSet&, const RefGaussianBasisSet&,
                 OneBodyIntIter* =0);
    ~OneBodyIntv2();
};

class OneBody3Intv2 : public OneBody3Int
{
  private:
    int same_center;
  protected:
    struct struct_centers* c1;
    struct struct_centers* c2;
  public:
    OneBody3Intv2(const RefGaussianBasisSet&, OneBodyIntIter* =0);
    OneBody3Intv2(const RefGaussianBasisSet&, const RefGaussianBasisSet&,
                  OneBodyIntIter* =0);
    ~OneBody3Intv2();
};

class GaussianOverlapIntv2 : public OneBodyIntv2
{
  public:
    GaussianOverlapIntv2(const RefGaussianBasisSet&, OneBodyIntIter* =0);
    GaussianOverlapIntv2(const RefGaussianBasisSet&,
                         const RefGaussianBasisSet&,
                         OneBodyIntIter* =0);
    ~GaussianOverlapIntv2();
    void compute_shell(int,int,double*);
};

class GaussianKineticIntv2 : public OneBodyIntv2
{
  public:
    GaussianKineticIntv2(const RefGaussianBasisSet&, OneBodyIntIter* =0);
    GaussianKineticIntv2(const RefGaussianBasisSet&,
                         const RefGaussianBasisSet&,
                         OneBodyIntIter* =0);
    ~GaussianKineticIntv2();
    void compute_shell(int,int,double*);
};

class GaussianPointChargeIntv2 : public OneBodyIntv2
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

class GaussianNuclearIntv2 : public GaussianPointChargeIntv2
{
  public:
    GaussianNuclearIntv2(const RefGaussianBasisSet&, OneBodyIntIter* =0);
    GaussianNuclearIntv2(PointBag_double *charges,
                         const RefGaussianBasisSet&,
                         const RefGaussianBasisSet&,
                         OneBodyIntIter* =0);

    ~GaussianNuclearIntv2();
};

class GaussianEfieldDotVectorIntv2: public OneBodyIntv2
{
  private:
    double *buffer3_; // a larger buffer is needed than that provided
    double position_[3];
    double vector_[3];
  public:
    GaussianEfieldDotVectorIntv2(const RefGaussianBasisSet&,
                                 double *postion = 0,
                                 double *vector = 0,
                                 OneBodyIntIter* =0);
    GaussianEfieldDotVectorIntv2(const RefGaussianBasisSet&,
                                 const RefGaussianBasisSet&,
                                 double *postion = 0,
                                 double *vector = 0,
                                 OneBodyIntIter* =0);
    ~GaussianEfieldDotVectorIntv2();
    void position(const double*);
    void vector(const double*);
    void compute_shell(int,int,double*);
};

class GaussianDipoleIntv2: public OneBody3Intv2
{
  private:
    double origin_[3];
  public:
    GaussianDipoleIntv2(const RefGaussianBasisSet&,
                        const double *origin = 0,
                        OneBodyIntIter* =0);
    GaussianDipoleIntv2(const RefGaussianBasisSet&,
                        const RefGaussianBasisSet&,
                        const double *origin = 0,
                        OneBodyIntIter* =0);
    ~GaussianDipoleIntv2();
    void origin(const double*);
    void compute_shell(int,int,double*);
};

////////////////////////////////////////////////////////////////////////////

class IntegralV2 : public Integral {
#   define CLASSNAME IntegralV2
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    IntegralV2();
    IntegralV2(StateIn&);
    IntegralV2(const RefKeyVal&);
    
    void save_data_state(StateOut&);
    
    CartesianIter * new_cartesian_iter(int);
    RedundantCartesianIter * new_redundant_cartesian_iter(int);
    RedundantCartesianSubIter * new_redundant_cartesian_sub_iter(int);
    SphericalTransformIter * new_spherical_transform_iter(int, int=0);
    
    RefSCElementOp overlap_op(const RefGaussianBasisSet&, OneBodyIntIter* = 0);
    RefSCElementOp overlap_op(const RefGaussianBasisSet&,
                              const RefGaussianBasisSet&, OneBodyIntIter* = 0);

    RefSCElementOp kinetic_op(const RefGaussianBasisSet&, OneBodyIntIter* = 0);
    RefSCElementOp kinetic_op(const RefGaussianBasisSet&,
                              const RefGaussianBasisSet&, OneBodyIntIter* = 0);

    RefSCElementOp point_charge_op(PointBag_double*,
                                   const RefGaussianBasisSet&,
                                   OneBodyIntIter* = 0);
    RefSCElementOp point_charge_op(PointBag_double*,
                                   const RefGaussianBasisSet&,
                                   const RefGaussianBasisSet&,
                                   OneBodyIntIter* = 0);

    RefSCElementOp nuclear_op(const RefGaussianBasisSet&, OneBodyIntIter* = 0);
    RefSCElementOp nuclear_op(PointBag_double*, const RefGaussianBasisSet&,
                              const RefGaussianBasisSet&, OneBodyIntIter* = 0);

    RefSCElementOp efield_dot_vector_op(const RefGaussianBasisSet&,
                                        double *position = 0,
                                        double *vector = 0,
                                        OneBodyIntIter* = 0);
    RefSCElementOp efield_dot_vector_op(const RefGaussianBasisSet&,
                                        const RefGaussianBasisSet&,
                                        double *position = 0,
                                        double *vector = 0,
                                        OneBodyIntIter* = 0);

    RefSCElementOp3 dipole_op(const RefGaussianBasisSet&, double *origin = 0,
                              OneBodyIntIter* = 0);
    RefSCElementOp3 dipole_op(const RefGaussianBasisSet&,
                              const RefGaussianBasisSet&,
                              double *origin = 0, OneBodyIntIter* = 0);
};

#endif
