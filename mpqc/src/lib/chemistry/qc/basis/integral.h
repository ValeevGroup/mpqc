
#ifndef _chemistry_qc_basis_integral_h
#define _chemistry_qc_basis_integral_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/state/state.h>

class SymmetryOperation;
class RefPetiteList;
class RefOneBodyInt;
class RefOneBodyIntIter;
class RefOneBodyDerivInt;
class RefTwoBodyInt;
class RefTwoBodyDerivInt;
class RefSymmSCMatrix;
class ShellRotation;
class CartesianIter;
class RedundantCartesianIter;
class RedundantCartesianSubIter;
class SphericalTransformIter;
class PointBag_double;

SavableState_REF_fwddec(SCElementOp)
SavableState_REF_fwddec(SCElementOp3)
SavableState_REF_fwddec(GaussianBasisSet)


// some useful things to have that depend on the underlying integrals
// package

class Integral : public SavableState {
#   define CLASSNAME Integral
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    Integral();

  public:
    Integral(StateIn&);
    Integral(const RefKeyVal&);
    
    void save_data_state(StateOut&);

    RefPetiteList petite_list(const RefGaussianBasisSet&);
    ShellRotation shell_rotation(int am, SymmetryOperation&, int pure=0);


    ///////////////////////////////////////////////////////////////////////
    // the following must be defined in the specific integral package

    virtual CartesianIter * new_cartesian_iter(int) =0;
    virtual RedundantCartesianIter * new_redundant_cartesian_iter(int) =0;
    virtual RedundantCartesianSubIter *
                                 new_redundant_cartesian_sub_iter(int) =0;
    virtual SphericalTransformIter *
                              new_spherical_transform_iter(int, int=0) =0;
    
    virtual RefOneBodyInt overlap_int(const RefGaussianBasisSet&) =0;
    virtual RefOneBodyInt overlap_int(const RefGaussianBasisSet&,
                                      const RefGaussianBasisSet&) =0;
    
    virtual RefOneBodyInt kinetic_int(const RefGaussianBasisSet&) =0;
    virtual RefOneBodyInt kinetic_int(const RefGaussianBasisSet&,
                                      const RefGaussianBasisSet&) =0;
    

    virtual RefOneBodyInt point_charge_int(PointBag_double*,
                                           const RefGaussianBasisSet&) =0;
    virtual RefOneBodyInt point_charge_int(PointBag_double*,
                                           const RefGaussianBasisSet&,
                                           const RefGaussianBasisSet&) =0;

    virtual RefOneBodyInt nuclear_int(const RefGaussianBasisSet&) =0;
    virtual RefOneBodyInt nuclear_int(PointBag_double*,
                                      const RefGaussianBasisSet&,
                                      const RefGaussianBasisSet&) =0;

    virtual RefOneBodyInt efield_dot_vector_int(const RefGaussianBasisSet&,
                                                double *position = 0,
                                                double *vector = 0) =0;
    virtual RefOneBodyInt efield_dot_vector_int(const RefGaussianBasisSet&,
                                                const RefGaussianBasisSet&,
                                                double *position = 0,
                                                double *vector = 0) =0;

    virtual RefOneBodyInt dipole_int(const RefGaussianBasisSet&,
                                     double *origin = 0) =0;
    virtual RefOneBodyInt dipole_int(const RefGaussianBasisSet&,
                                     const RefGaussianBasisSet&,
                                     double *origin = 0) =0;

    virtual RefOneBodyDerivInt deriv_int(const RefGaussianBasisSet&) =0;
    virtual RefOneBodyDerivInt deriv_int(const RefGaussianBasisSet&,
                                         const RefGaussianBasisSet&) =0;
    
    virtual RefTwoBodyInt two_body_int(const RefGaussianBasisSet&) =0;
    virtual RefTwoBodyInt two_body_int(const RefGaussianBasisSet&,
                                       const RefGaussianBasisSet&,
                                       const RefGaussianBasisSet&,
                                       const RefGaussianBasisSet&) =0;
    
    virtual RefTwoBodyDerivInt two_body_deriv_int(const RefGaussianBasisSet&) =0;
    virtual RefTwoBodyDerivInt two_body_deriv_int(const RefGaussianBasisSet&,
                                                  const RefGaussianBasisSet&,
                                                  const RefGaussianBasisSet&,
                                                  const RefGaussianBasisSet&) =0;
};
SavableState_REF_dec(Integral);

#endif
