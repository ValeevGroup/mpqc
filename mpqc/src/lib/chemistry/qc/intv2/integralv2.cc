
#include <stdio.h>

#include <chemistry/qc/intv2/integralv2.h>
#include <chemistry/qc/intv2/cartiterv2.h>
#include <chemistry/qc/intv2/tformv2.h>
#include <chemistry/qc/intv2/obintv2.h>
#include <chemistry/qc/intv2/tbintv2.h>

#define CLASSNAME IntegralV2
#define PARENTS public Integral
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
IntegralV2::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = Integral::_castdown(cd);
  return do_castdowns(casts,cd);
}

IntegralV2::IntegralV2()
{
}

IntegralV2::IntegralV2(StateIn& s) :
  Integral(s)
{
}

IntegralV2::IntegralV2(const RefKeyVal& k) :
  Integral(k)
{
}

void
IntegralV2::save_data_state(StateOut& s)
{
  Integral::save_data_state(s);
}

CartesianIter *
IntegralV2::new_cartesian_iter(int l)
{
  return new CartesianIterV2(l);
}

RedundantCartesianIter *
IntegralV2::new_redundant_cartesian_iter(int l)
{
  return new RedundantCartesianIterV2(l);
}

RedundantCartesianSubIter *
IntegralV2::new_redundant_cartesian_sub_iter(int l)
{
  return new RedundantCartesianSubIterV2(l);
}

SphericalTransformIter *
IntegralV2::new_spherical_transform_iter(int l, int inv)
{
  return new SphericalTransformIterV2(l, inv);
}

RefOneBodyInt
IntegralV2::overlap_int(const RefGaussianBasisSet& gbs)
{
  return new GaussianOverlapIntv2(gbs);
}

RefOneBodyInt
IntegralV2::overlap_int(const RefGaussianBasisSet& gbs1,
                        const RefGaussianBasisSet& gbs2)
{
  return new GaussianOverlapIntv2(gbs1, gbs2);
}

RefOneBodyInt
IntegralV2::kinetic_int(const RefGaussianBasisSet& gbs)
{
  return new GaussianKineticIntv2(gbs);
}

RefOneBodyInt
IntegralV2::kinetic_int(const RefGaussianBasisSet& gbs1,
                        const RefGaussianBasisSet& gbs2)
{
  return new GaussianKineticIntv2(gbs1, gbs2);
}

RefOneBodyInt
IntegralV2::point_charge_int(PointBag_double* pb,
                             const RefGaussianBasisSet& gbs)
{
  return new GaussianPointChargeIntv2(pb, gbs);
}

RefOneBodyInt
IntegralV2::point_charge_int(PointBag_double *pb,
                             const RefGaussianBasisSet& gbs1,
                             const RefGaussianBasisSet& gbs2)
{
  return new GaussianPointChargeIntv2(pb, gbs1, gbs2);
}

RefOneBodyInt
IntegralV2::nuclear_int(const RefGaussianBasisSet& gbs)
{
  return new GaussianNuclearIntv2(gbs);
}

RefOneBodyInt
IntegralV2::nuclear_int(PointBag_double *pb,
                        const RefGaussianBasisSet& gbs1,
                        const RefGaussianBasisSet& gbs2)
{
  return new GaussianNuclearIntv2(pb, gbs1, gbs2);
}

RefOneBodyInt
IntegralV2::efield_dot_vector_int(const RefGaussianBasisSet& gbs,
                                 double *position,
                                 double *vector)
{
  return new GaussianEfieldDotVectorIntv2(gbs, position, vector);
}

RefOneBodyInt
IntegralV2::efield_dot_vector_int(const RefGaussianBasisSet& gbs1,
                                 const RefGaussianBasisSet& gbs2,
                                 double *position,
                                 double *vector)
{
  return new GaussianEfieldDotVectorIntv2(gbs1, gbs2, position, vector);
}

RefOneBodyInt
IntegralV2::dipole_int(const RefGaussianBasisSet& gbs,
                      double *origin)
{
  return new GaussianDipoleIntv2(gbs, origin);
}

RefOneBodyInt
IntegralV2::dipole_int(const RefGaussianBasisSet& gbs1,
                       const RefGaussianBasisSet& gbs2,
                       double *origin)
{
  return new GaussianDipoleIntv2(gbs1, gbs2, origin);
}

RefOneBodyDerivInt
IntegralV2::deriv_int(const RefGaussianBasisSet& gbs)
{
  return new OneBodyDerivIntv2(gbs);
}

RefOneBodyDerivInt
IntegralV2::deriv_int(const RefGaussianBasisSet& gbs1,
                      const RefGaussianBasisSet& gbs2)
{
  return new OneBodyDerivIntv2(gbs1, gbs2);
}

RefTwoBodyInt
IntegralV2::two_body_int(const RefGaussianBasisSet& gbs)
{
  return new TwoBodyIntV2(gbs);
}

RefTwoBodyInt
IntegralV2::two_body_int(const RefGaussianBasisSet& gbs1,
                         const RefGaussianBasisSet& gbs2,
                         const RefGaussianBasisSet& gbs3,
                         const RefGaussianBasisSet& gbs4)
{
  return new TwoBodyIntV2(gbs1, gbs2, gbs3, gbs4);
}

RefTwoBodyDerivInt
IntegralV2::two_body_deriv_int(const RefGaussianBasisSet& gbs)
{
  return new TwoBodyDerivIntV2(gbs);
}

RefTwoBodyDerivInt
IntegralV2::two_body_deriv_int(const RefGaussianBasisSet& gbs1,
                         const RefGaussianBasisSet& gbs2,
                         const RefGaussianBasisSet& gbs3,
                         const RefGaussianBasisSet& gbs4)
{
  return new TwoBodyDerivIntV2(gbs1, gbs2, gbs3, gbs4);
}
