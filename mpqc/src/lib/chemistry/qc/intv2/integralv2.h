
// these provide integrals using the libintv2 routines

#ifndef _chemistry_qc_intv2_integralv2_h
#define _chemistry_qc_intv2_integralv2_h

#include <math/topology/pointbag.h>

#include <chemistry/qc/basis/integral.h>

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
    
    RefOneBodyInt overlap_int(const RefGaussianBasisSet&);
    RefOneBodyInt overlap_int(const RefGaussianBasisSet&,
                              const RefGaussianBasisSet&);

    RefOneBodyInt kinetic_int(const RefGaussianBasisSet&);
    RefOneBodyInt kinetic_int(const RefGaussianBasisSet&,
                              const RefGaussianBasisSet&);

    RefOneBodyInt point_charge_int(PointBag_double*,
                                   const RefGaussianBasisSet&);
    RefOneBodyInt point_charge_int(PointBag_double*,
                                   const RefGaussianBasisSet&,
                                   const RefGaussianBasisSet&);

    RefOneBodyInt nuclear_int(const RefGaussianBasisSet&);
    RefOneBodyInt nuclear_int(PointBag_double*, const RefGaussianBasisSet&,
                              const RefGaussianBasisSet&);

    RefOneBodyInt efield_dot_vector_int(const RefGaussianBasisSet&,
                                        double *position = 0,
                                        double *vector = 0);
    RefOneBodyInt efield_dot_vector_int(const RefGaussianBasisSet&,
                                        const RefGaussianBasisSet&,
                                        double *position = 0,
                                        double *vector = 0);

    RefOneBodyInt dipole_int(const RefGaussianBasisSet&, double *origin = 0);
    RefOneBodyInt dipole_int(const RefGaussianBasisSet&,
                             const RefGaussianBasisSet&,
                             double *origin =0);

    RefOneBodyDerivInt deriv_int(const RefGaussianBasisSet&);
    RefOneBodyDerivInt deriv_int(const RefGaussianBasisSet&,
                                 const RefGaussianBasisSet&);

    RefTwoBodyInt two_body_int(const RefGaussianBasisSet&);
    RefTwoBodyInt two_body_int(const RefGaussianBasisSet&,
                               const RefGaussianBasisSet&,
                               const RefGaussianBasisSet&,
                               const RefGaussianBasisSet&);

    RefTwoBodyDerivInt two_body_deriv_int(const RefGaussianBasisSet&);
    RefTwoBodyDerivInt two_body_deriv_int(const RefGaussianBasisSet&,
                                          const RefGaussianBasisSet&,
                                          const RefGaussianBasisSet&,
                                          const RefGaussianBasisSet&);
};


#endif
