
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
IntegralV2::overlap()
{
  return new GaussianOverlapIntv2(bs1_, bs2_);
}

RefOneBodyInt
IntegralV2::kinetic()
{
  return new GaussianKineticIntv2(bs1_, bs2_);
}

RefOneBodyInt
IntegralV2::point_charge(const RefPointChargeData& dat)
{
  return new GaussianPointChargeIntv2(bs1_, bs2_, dat);
}

RefOneBodyInt
IntegralV2::nuclear()
{
  return new GaussianNuclearIntv2(bs1_, bs2_);
}

RefOneBodyInt
IntegralV2::efield_dot_vector(const RefEfieldDotVectorData&dat)
{
  return new GaussianEfieldDotVectorIntv2(bs1_, bs2_, dat);
}

RefOneBodyInt
IntegralV2::dipole(const RefDipoleData& dat)
{
  return new GaussianDipoleIntv2(bs1_, bs2_, dat);
}

RefOneBodyDerivInt
IntegralV2::deriv()
{
  return new OneBodyDerivIntv2(bs1_, bs2_);
}

RefTwoBodyInt
IntegralV2::electron_repulsion()
{
  return new TwoBodyIntV2(bs1_, bs2_, bs3_, bs4_);
}

RefTwoBodyDerivInt
IntegralV2::electron_repulsion_deriv()
{
  return new TwoBodyDerivIntV2(bs1_, bs2_, bs3_, bs4_);
}
