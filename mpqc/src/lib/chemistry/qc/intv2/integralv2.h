
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
    
    RefOneBodyInt overlap();

    RefOneBodyInt kinetic();

    RefOneBodyInt point_charge(const RefPointChargeData& =0);

    RefOneBodyInt nuclear();

    RefOneBodyInt efield_dot_vector(const RefEfieldDotVectorData& =0);

    RefOneBodyInt dipole(const RefDipoleData& =0);

    RefOneBodyDerivInt deriv();

    RefTwoBodyInt electron_repulsion();

    RefTwoBodyDerivInt electron_repulsion_deriv();
};


#endif
