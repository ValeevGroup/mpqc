//
// intv3.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <util/state/stateio.h>
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
  initialize_transforms();
}

IntegralV3::IntegralV3(StateIn& s) :
  Integral(s)
{
  initialize_transforms();
}

IntegralV3::IntegralV3(const RefKeyVal& k) :
  Integral(k)
{
  initialize_transforms();
}

void
IntegralV3::save_data_state(StateOut& s)
{
  Integral::save_data_state(s);
}

IntegralV3::~IntegralV3()
{
  free_transforms();
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
IntegralV3::new_spherical_transform_iter(int l, int inv, int subl)
{
  if (l>maxl_ || l<0) {
      cerr << "IntegralV3::new_spherical_transform_iter: bad l" << endl;
      abort();
    }
  if (subl == -1) subl = l;
  if (subl < 0 || subl > l || (l-subl)%2 != 0) {
      cerr << "IntegralV3::new_spherical_transform_iter: bad subl" << endl;
      abort();
    }
  if (inv) {
      return new SphericalTransformIter(ist_[l][(l-subl)/2]);
    }
  return new SphericalTransformIter(st_[l][(l-subl)/2]);
}

const SphericalTransform *
IntegralV3::spherical_transform(int l, int inv, int subl)
{
  if (l>maxl_ || l<0) {
      cerr << "IntegralV3::spherical_transform_iter: bad l" << endl;
      abort();
    }
  if (subl == -1) subl = l;
  if (subl < 0 || subl > l || (l-subl)%2 != 0) {
      cerr << "IntegralV3::spherical_transform_iter: bad subl" << endl;
      abort();
    }
  if (inv) {
      return ist_[l][(l-subl)/2];
    }
  return st_[l][(l-subl)/2];
}

RefOneBodyInt
IntegralV3::overlap()
{
  return new OneBodyIntV3(this, bs1_, bs2_, &Int1eV3::overlap);
}

RefOneBodyInt
IntegralV3::kinetic()
{
  return new OneBodyIntV3(this, bs1_, bs2_, &Int1eV3::kinetic);
}

RefOneBodyInt
IntegralV3::nuclear()
{
  return new OneBodyIntV3(this, bs1_, bs2_, &Int1eV3::nuclear);
}

RefOneBodyInt
IntegralV3::hcore()
{
  return new OneBodyIntV3(this, bs1_, bs2_, &Int1eV3::hcore);
}

RefOneBodyInt
IntegralV3::point_charge(const RefPointChargeData& dat)
{
  return new PointChargeIntV3(this, bs1_, bs2_, dat);
}

RefOneBodyInt
IntegralV3::efield_dot_vector(const RefEfieldDotVectorData&dat)
{
  return new EfieldDotVectorIntV3(this, bs1_, bs2_, dat);
}

RefOneBodyInt
IntegralV3::dipole(const RefDipoleData& dat)
{
  return new DipoleIntV3(this, bs1_, bs2_, dat);
}

RefOneBodyDerivInt
IntegralV3::overlap_deriv()
{
  return new OneBodyDerivIntV3(this, bs1_, bs2_, &Int1eV3::overlap_1der);
}

RefOneBodyDerivInt
IntegralV3::kinetic_deriv()
{
  return new OneBodyDerivIntV3(this, bs1_, bs2_, &Int1eV3::kinetic_1der);
}

RefOneBodyDerivInt
IntegralV3::nuclear_deriv()
{
  return new OneBodyDerivIntV3(this, bs1_, bs2_, &Int1eV3::nuclear_1der);
}

RefOneBodyDerivInt
IntegralV3::hcore_deriv()
{
  return new OneBodyDerivIntV3(this, bs1_, bs2_, &Int1eV3::hcore_1der);
}

RefTwoBodyInt
IntegralV3::electron_repulsion()
{
  return new TwoBodyIntV3(this, bs1_, bs2_, bs3_, bs4_, storage_);
}

RefTwoBodyDerivInt
IntegralV3::electron_repulsion_deriv()
{
  return new TwoBodyDerivIntV3(this, bs1_, bs2_, bs3_, bs4_, storage_);
}

void
IntegralV3::set_basis(const RefGaussianBasisSet &b1,
                      const RefGaussianBasisSet &b2,
                      const RefGaussianBasisSet &b3,
                      const RefGaussianBasisSet &b4)
{
  free_transforms();
  Integral::set_basis(b1,b2,b3,b4);
  initialize_transforms();
}

void
IntegralV3::free_transforms()
{
  int i,j;
  for (i=0; i<=maxl_; i++) {
      for (j=0; j<=i/2; j++) {
          delete st_[i][j];
          delete ist_[i][j];
        }
      delete[] st_[i];
      delete[] ist_[i];
    }
  delete[] st_;
  delete[] ist_;
  st_ = 0;
  ist_ = 0;
}

void
IntegralV3::initialize_transforms()
{
  maxl_ = -1;
  int maxam;
  maxam = bs1_.nonnull()?bs1_->max_angular_momentum():-1;
  if (maxl_ < maxam) maxl_ = maxam;
  maxam = bs2_.nonnull()?bs2_->max_angular_momentum():-1;
  if (maxl_ < maxam) maxl_ = maxam;
  maxam = bs3_.nonnull()?bs3_->max_angular_momentum():-1;
  if (maxl_ < maxam) maxl_ = maxam;
  maxam = bs4_.nonnull()?bs4_->max_angular_momentum():-1;
  if (maxl_ < maxam) maxl_ = maxam;

  st_ = new SphericalTransformV3**[maxl_+1];
  ist_ = new ISphericalTransformV3**[maxl_+1];;
  int i,j;
  for (i=0; i<=maxl_; i++) {
      st_[i] = new SphericalTransformV3*[i/2+1];
      ist_[i] = new ISphericalTransformV3*[i/2+1];
      for (j=0; j<=i/2; j++) {
          st_[i][j] = new SphericalTransformV3(i,i-2*j);
          ist_[i][j] = new ISphericalTransformV3(i,i-2*j);
        }
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
