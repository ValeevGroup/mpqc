
#include <chemistry/qc/intv3/intv3.h>
#include <chemistry/qc/intv3/cartitv3.h>
#include <chemistry/qc/intv3/tformv3.h>
#include <chemistry/qc/intv3/obintv3.h>
#include <chemistry/qc/intv3/tbintv3.h>

#define CLASSNAME IntegralV3
#define PARENTS public Integral
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
IntegralV3::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = Integral::_castdown(cd);
  return do_castdowns(casts,cd);
}

IntegralV3::IntegralV3(const RefGaussianBasisSet &b1,
                       const RefGaussianBasisSet &b2,
                       const RefGaussianBasisSet &b3,
                       const RefGaussianBasisSet &b4):
  Integral(b1,b2,b3,b4)
{
}

IntegralV3::IntegralV3(StateIn& s) :
  Integral(s)
{
}

IntegralV3::IntegralV3(const RefKeyVal& k) :
  Integral(k)
{
}

void
IntegralV3::save_data_state(StateOut& s)
{
  Integral::save_data_state(s);
}

CartesianIter *
IntegralV3::new_cartesian_iter(int l)
{
  return new CartesianIterV3(l);
}

RedundantCartesianIter *
IntegralV3::new_redundant_cartesian_iter(int l)
{
  return new RedundantCartesianIterV3(l);
}

RedundantCartesianSubIter *
IntegralV3::new_redundant_cartesian_sub_iter(int l)
{
  return new RedundantCartesianSubIterV3(l);
}

SphericalTransformIter *
IntegralV3::new_spherical_transform_iter(int l, int inv)
{
  return new SphericalTransformIterV3(l, inv);
}

RefOneBodyInt
IntegralV3::overlap()
{
  return new OneBodyIntV3(bs1_, bs2_, Int1eV3::overlap);
}

RefOneBodyInt
IntegralV3::kinetic()
{
  return new OneBodyIntV3(bs1_, bs2_, Int1eV3::kinetic);
}

RefOneBodyInt
IntegralV3::nuclear()
{
  return new OneBodyIntV3(bs1_, bs2_, Int1eV3::nuclear);
}

RefOneBodyInt
IntegralV3::hcore()
{
  return new OneBodyIntV3(bs1_, bs2_, Int1eV3::hcore);
}

RefOneBodyInt
IntegralV3::point_charge(const RefPointChargeData& dat)
{
  return new PointChargeIntV3(bs1_, bs2_, dat);
}

RefOneBodyInt
IntegralV3::efield_dot_vector(const RefEfieldDotVectorData&dat)
{
  return new EfieldDotVectorIntV3(bs1_, bs2_, dat);
}

RefOneBodyInt
IntegralV3::dipole(const RefDipoleData& dat)
{
  return new DipoleIntV3(bs1_, bs2_, dat);
}

RefOneBodyDerivInt
IntegralV3::overlap_deriv()
{
  return new OneBodyDerivIntV3(bs1_, bs2_, Int1eV3::overlap_1der);
}

RefOneBodyDerivInt
IntegralV3::kinetic_deriv()
{
  return new OneBodyDerivIntV3(bs1_, bs2_, Int1eV3::kinetic_1der);
}

RefOneBodyDerivInt
IntegralV3::nuclear_deriv()
{
  return new OneBodyDerivIntV3(bs1_, bs2_, Int1eV3::nuclear_1der);
}

RefOneBodyDerivInt
IntegralV3::hcore_deriv()
{
  return new OneBodyDerivIntV3(bs1_, bs2_, Int1eV3::hcore_1der);
}

RefTwoBodyInt
IntegralV3::electron_repulsion()
{
  return new TwoBodyIntV3(bs1_, bs2_, bs3_, bs4_, storage_);
}

RefTwoBodyDerivInt
IntegralV3::electron_repulsion_deriv()
{
  return new TwoBodyDerivIntV3(bs1_, bs2_, bs3_, bs4_, storage_);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
